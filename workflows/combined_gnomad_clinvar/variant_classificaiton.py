#variant_classification
import re
from enum import Enum
from typing import ClassVar, Dict, Optional, Tuple
"""Module for HGVS tokenization."""
import re
from typing import Optional

from cool_seq_tool.schemas import AnnotationLayer

from variation.schemas.token_response_schema import HgvsToken
from variation.tokenizers.tokenizer import Tokenizer

coord_types: ClassVar[Dict[str, str]] = {
        k: v.value for k, v in AnnotationLayer.__members__.items()
    }
class AnnotationLayer(str, Enum):
    """Create enum for supported annotation layers"""

    PROTEIN = "p"
    CDNA = "c"
    GENOMIC = "g"
#idk

class HGVS(Tokenizer):
    """The HGVS tokenizer class."""

    splitter = re.compile(
        r"^(?P<accession>(NC_|NM_|NP_|ENSP|ENST)[^:\s]+):(?P<coordinate>[cgnpr])\.(?P<change>\S+)$"
    )

    def match(self, input_string: str) -> Optional[HgvsToken]:
        """Return HGVS token matches from input string.

        :param input_string: The input string to match
        :return: `HgvsToken` if HGVS match was found, else `None`
        """
        match = self.splitter.match(input_string)
        if match:
            match_dict = match.groupdict()

            return HgvsToken(
                token=input_string,
                input_string=input_string,
                accession=match_dict["accession"],
                coordinate_type=AnnotationLayer(match_dict["coordinate"]),
                change=match_dict["change"],
            )

        return None
    
class GnomadVCF(Tokenizer):
    """The gnomad VCF tokenizer class"""

    splitter = re.compile(
        r"^(chr|chromosome)?(?P<chromosome>([1-9]|[1][0-9]|[2][0-2]|X|Y))-"
        r"(?P<pos>[1-9]\d*)-(?P<ref>[actg]+)-(?P<alt>[actg]+)$",
        re.IGNORECASE,
    )

    def match(self, input_string: str):
        """Return a GnomadVCFToken if a match exists.

        :param input_string: The input string to match
        :return: `Token` if gnomAD VCF match was found, else `None`
        """
        match = self.splitter.match(input_string)
        if match:
            match_dict = match.groupdict()
            chromosome = match_dict["chromosome"].upper()
            pos = int(match_dict["pos"])
            ref = match_dict["ref"].upper()
            alt = match_dict["alt"].upper()

            return GnomadVcfToken(
                token=f"{chromosome}-{pos}-{ref}-{alt}",
                input_string=input_string,
                chromosome=chromosome,
                pos=pos,
                ref=ref,
                alt=alt,
            )

        return None



"""Module for Tokenization."""
from abc import ABC, abstractmethod


import AnnotationLayer

from variation.schemas.token_response_schema import Token




"""A module for tokenization."""
from typing import List

from variation.schemas.token_response_schema import Token, TokenType
from variation.tokenizers import (
    HGVS,
    GnomadVCF,

)
from variation.tokenizers.tokenizer import Tokenizer



class Tokenize:
    """The tokenize class."""

    def __init__(self, gene_symbol: GeneSymbol) -> None:
        """Initialize the tokenize class."""
        self.tokenizers: List[Tokenizer] = [
            HGVS(),
            GnomadVCF(),
            #FreeTextCategorical(),
            # Substitution
            ProteinSubstitution(),
            GenomicSubstitution(),
            CdnaSubstitution(),
            # Reference Agree
            ProteinReferenceAgree(),
            CdnaGenomicReferenceAgree(),
            # Delins
            ProteinDelIns(),
            CdnaDelIns(),
            GenomicDelIns(),
            # Deletion
            ProteinDeletion(),
            CdnaDeletion(),
            GenomicDeletion(),
            # Insertion
            ProteinInsertion(),
            CdnaInsertion(),
            GenomicInsertion(),
            # Duplication
            GenomicDuplication(),
        ]

    def perform(self, search_string: str, warnings: List[str]) -> List[Token]:
        """Return a list of tokens for a given search string

        :param search_string: The input string to search on
        :param warnings: List of warnings
        :return: A list of tokens found
        """
        terms = search_string.split()

        tokens: List[Token] = []
        for term in terms:
            if not term:
                continue

            matched = False
            for tokenizer in self.tokenizers:
                res = tokenizer.match(term)
                if res:
                    if isinstance(res, List):
                        for r in res:
                            tokens.append(r)
                            if not matched:
                                matched = True
                    else:
                        tokens.append(res)
                        matched = True
                        break

            if not matched:
                warnings.append(f"Unable to tokenize: {term}")
                tokens.append(
                    Token(token=term, token_type=TokenType.UNKNOWN, input_string=term)
                )

        return tokens
    

class Tokenizer(ABC):
    """The tokenizer class."""

    coord_types: ClassVar[Dict[str, str]] = {
        k: v.value for k, v in AnnotationLayer.__members__.items()
    }

    @abstractmethod
    def match(self, input_string: str) -> Optional[Token]:
        """Return tokens that match the input string.

        :param input_string: Input string
        :return: Token if match was found
        """
        raise NotImplementedError

    def strip_coord_prefix(
        self, input_string: str, match_coord_type: Optional[AnnotationLayer] = None
    ) -> Tuple[Optional[AnnotationLayer], Optional[str]]:
        """Strip parentheses and coordinate type from string

        :param input_string: Input string
        :param match_coord_type: If set, the input string must have the prefix
            corresponding to this value to succeed. If this is not set, will attempt
            to find the first match of a prefix and use that as the coordinate type.
        :return: Tuple containing coordinate type for input string and stripped string,
            if successful.
        """
        coord_type = None
        stripped_str = None

        def _strip(
            coord_type: str,
            string: str,
            match_coord_type: Optional[AnnotationLayer] = None,
        ) -> str:
            """Strip parentheses and coordinate type from string

            :param input_string: Input string
            :param match_coord_type: If set, the input string must have the prefix
                corresponding to this value to succeed
            :return: Stripped string
            """
            if string.startswith(
                (f"({coord_type}.", f"{coord_type}.(")
            ) and string.endswith(")"):
                string = string[3:-1]
            elif string.startswith(f"{coord_type}."):
                string = string[2:]
            elif string[0] == "(" and string[-1] == ")":
                string = string[1:-1]
            else:
                if match_coord_type:
                    string = None

            return string

        if match_coord_type:
            coord_type = match_coord_type
            stripped_str = _strip(coord_type.value, input_string, match_coord_type)
        else:
            for k, v in self.coord_types.items():
                if f"{v}." in input_string:
                    coord_type = AnnotationLayer[k]
                    stripped_str = _strip(v, input_string)
                    break

        return coord_type, stripped_str
    


    


    
"""Module for classification."""
from typing import ClassVar, List, Optional

from variation.classifiers import (
    GnomadVcfClassifier,
    HgvsClassifier,
)
from variation.classifiers.classifier import Classifier
from variation.schemas.classification_response_schema import Classification
from variation.schemas.token_response_schema import Token, TokenType


class Classify:
    """The classify class."""

    hgvs_classifier = HgvsClassifier()
    gnomad_vcf_classifier = GnomadVcfClassifier()

    def perform(self, tokens: List[Token]) -> Optional[Classification]:
        """Classify a list of tokens.

        :param tokens: List of tokens found
        :return: Classification for a list of tokens if found
        """
        classification = None

        if len(tokens) == 1:
            token_type = tokens[0].token_type

            if token_type == TokenType.HGVS:
                classification = self.hgvs_classifier.match(tokens[0])
            elif token_type == TokenType.GNOMAD_VCF:
                classification = self.gnomad_vcf_classifier.match(tokens[0])
        else:
            raise TypeError

        return classification
    

"""A module for the gnomAD VCF Classifier"""
from typing import List, Optional, Union

from variation.classifiers.classifier import Classifier
from variation.schemas.classification_response_schema import (
    GenomicDeletionClassification,
    GenomicDelInsClassification,
    GenomicInsertionClassification,
    GenomicReferenceAgreeClassification,
    GenomicSubstitutionClassification,
    Nomenclature,
)
from variation.schemas.token_response_schema import GnomadVcfToken, TokenType


class GnomadVcfClassifier(Classifier):
    """The gnomAD VCF Classifier"""

    def exact_match_candidates(self) -> List[List[TokenType]]:
        """Return the token match candidates for the gnomad vcf classification.

        :return: List of list of tokens, where order matters, that represent a gnomad
        vcf classification.
        """
        return [[TokenType.GNOMAD_VCF]]

    def match(
        self, token: GnomadVcfToken
    ) -> Optional[
        Union[
            GenomicReferenceAgreeClassification,
            GenomicSubstitutionClassification,
            GenomicInsertionClassification,
            GenomicDeletionClassification,
        ]
    ]:
        """Return the genomic classification (either reference agree, substitution,
        insertion, or deletion) from a gnomad vcf token.
        Currently only support simple genomic variation.

        :param token: gnomad vcf token
        :return: The corresponding genomic classification for the gnomad vcf token if
            simple variation change. Else, return `None`
        """
        params = {"matching_tokens": [token], "nomenclature": Nomenclature.GNOMAD_VCF}

        ref = token.ref
        alt = token.alt

        len_ref = len(ref)
        len_alt = len(alt)

        if len_ref == len_alt:
            # substitution
            params["pos"] = token.pos

            if ref == alt:
                return GenomicReferenceAgreeClassification(**params)

            params["ref"] = ref
            params["alt"] = alt

            return GenomicSubstitutionClassification(**params)

        # delins
        params["pos0"] = token.pos
        params["pos1"] = (params["pos0"] + len_ref) - 1
        if params["pos0"] == params["pos1"]:
            del params["pos1"]

        params["inserted_sequence"] = alt
        return GenomicDelInsClassification(**params)


"""Module for Classification methods."""
from abc import ABC, abstractmethod
from typing import List, Optional

from variation.schemas.classification_response_schema import Classification
from variation.schemas.token_response_schema import Token, TokenType


class Classifier(ABC):
    """The Classifier class."""

    @abstractmethod
    def match(self, tokens: List[Token]) -> Optional[Classification]:
        """Return the classification from a list of token matches.

        :param tokens: List of ordered tokens that are exact match candidates for a
            given classification
        :return: A classification for the list of matched tokens
        """

    @abstractmethod
    def exact_match_candidates(self) -> List[List[TokenType]]:
        """Return the token match candidates for a given classification.

        :return: List of list of tokens, where order matters, that represent a given
            classification.
        """

    def can_classify(self, tokens: List[Token]) -> bool:
        """Return whether or not a list of tokens can be classified by a given
        classification

        :param tokens: List of tokens found in an input query
        :return: `True` if a list of tokens matches the tokens needed, where order
            matters, to represent a given classification. `False`, otherwise.
        """
        token_types = [t.token_type for t in tokens]
        exact_matches: List[List[str]] = []

        for candidate in self.exact_match_candidates():
            if token_types == candidate:
                exact_matches.append(candidate)

        return len(exact_matches) == 1
    
"""A module for the HGVS Classifier."""
from re import Match, Pattern
from typing import Dict, List, Optional

from cool_seq_tool.schemas import AnnotationLayer

from variation.classifiers.classifier import Classifier
from variation.regex import (
    CDNA_REGEXPRS,
    GENOMIC_DEL_AMBIGUOUS_REGEXPRS,
    GENOMIC_DUP_AMBIGUOUS_REGEXPRS,
    GENOMIC_REGEXPRS,
    PROTEIN_REGEXPRS,
)
from variation.schemas.app_schemas import AmbiguousRegexType
from variation.schemas.classification_response_schema import (
    CdnaDeletionClassification,
    CdnaDelInsClassification,
    CdnaInsertionClassification,
    CdnaReferenceAgreeClassification,
    CdnaSubstitutionClassification,
    Classification,
    ClassificationType,
    GenomicDeletionAmbiguousClassification,
    GenomicDeletionClassification,
    GenomicDelInsClassification,
    GenomicDuplicationAmbiguousClassification,
    GenomicDuplicationClassification,
    GenomicInsertionClassification,
    GenomicReferenceAgreeClassification,
    GenomicSubstitutionClassification,
    Nomenclature,
    ProteinDeletionClassification,
    ProteinDelInsClassification,
    ProteinInsertionClassification,
    ProteinReferenceAgreeClassification,
    ProteinStopGainClassification,
    ProteinSubstitutionClassification,
)
from variation.schemas.token_response_schema import HgvsToken, TokenType
from variation.utils import get_ambiguous_type


class HgvsClassifier(Classifier):
    """The HGVS Classifier."""

    def exact_match_candidates(self) -> List[List[TokenType]]:
        """Return the token match candidates for the hgvs classification.

        :return: List of list of tokens, where order matters, that represent a hgvs
        classification.
        """
        return [[TokenType.HGVS]]

    def match(self, token: HgvsToken) -> Optional[Classification]:
        """Return the classification from a hgvs token using regex matches to determine
        the type of classification.

        :param token: hgvs token
        :return: The corresponding classification for the hgvs token if a regex match
            is found. Else, return `None`
        """
        classification = None
        params = {
            "matching_tokens": [token],
            "nomenclature": Nomenclature.HGVS,
            "ac": token.accession,
        }

        if token.coordinate_type == AnnotationLayer.GENOMIC:
            classification = self._genomic_classification(token, params)
            if not classification:
                # Try ambiguous
                classification = self._genomic_ambiguous_classification(token, params)
        elif token.coordinate_type == AnnotationLayer.CDNA:
            classification = self._cdna_classification(token, params)
        elif token.coordinate_type == AnnotationLayer.PROTEIN:
            classification = self._protein_classification(token, params)

        return classification

    @staticmethod
    def _regex_match(change: str, regex: Pattern) -> Optional[Match]:
        """Strip parentheses from `change` and return whether or not `change` matches
        the `regex`

        :param change: The alteration part of the hgvs expression
        :param regex: The pattern to match against
        :return: A regex match if found against pattern, else `None`
        """
        if change[0] == "(" and change[-1] == ")":
            match = regex.match(change[1:-1])
        else:
            match = regex.match(change)
        return match

    def _protein_classification(
        self, token: HgvsToken, params: Dict
    ) -> Optional[Classification]:
        """Determine if hgvs token matches regex checks and return corresponding protein
        classification if a match is found

        :param token: hgvs token
        :param params: Base fields for a classification. This will get mutated if a
            match is found.
        :return: Protein classification if hgvs token matches regex checks. Else, `None`
        """
        for regex, _, classification_type in PROTEIN_REGEXPRS:
            match = self._regex_match(token.change, regex)

            if match:
                match_dict = match.groupdict()
                params.update(match_dict)

                if classification_type == ClassificationType.PROTEIN_SUBSTITUTION:
                    params["pos"] = int(params["pos"])
                    if params["alt"] in {"Ter", "*"}:
                        params["alt"] = "*"
                        return ProteinStopGainClassification(**params)

                    return ProteinSubstitutionClassification(**params)

                if classification_type == ClassificationType.PROTEIN_REFERENCE_AGREE:
                    params["pos"] = int(params["pos"])
                    return ProteinReferenceAgreeClassification(**params)

                if classification_type == ClassificationType.PROTEIN_DELINS:
                    params["pos0"] = int(params["pos0"])
                    params["pos1"] = (
                        int(params["pos1"])
                        if params["pos1"] is not None
                        else params["pos1"]
                    )
                    return ProteinDelInsClassification(**params)

                if classification_type == ClassificationType.PROTEIN_DELETION:
                    params["pos0"] = int(params["pos0"])
                    params["pos1"] = (
                        int(params["pos1"])
                        if params["pos1"] is not None
                        else params["pos1"]
                    )
                    return ProteinDeletionClassification(**params)

                if classification_type == ClassificationType.PROTEIN_INSERTION:
                    params["pos0"] = int(params["pos0"])
                    params["pos1"] = (
                        int(params["pos1"])
                        if params["pos1"] is not None
                        else params["pos1"]
                    )
                    return ProteinInsertionClassification(**params)

        return None

    def _cdna_classification(
        self, token: HgvsToken, params: Dict
    ) -> Optional[Classification]:
        """Determine if hgvs token matches regex checks and return corresponding cdna
        classification if a match is found

        :param token: hgvs token
        :param params: Base fields for a classification. This will get mutated if a
            match is found.
        :return: cdna classification if hgvs token matches regex checks. Else, `None`
        """
        for regex, _, classification_type in CDNA_REGEXPRS:
            match = self._regex_match(token.change, regex)

            if match:
                match_dict = match.groupdict()
                params.update(match_dict)

                if classification_type == ClassificationType.CDNA_SUBSTITUTION:
                    params["pos"] = int(params["pos"])
                    return CdnaSubstitutionClassification(**params)

                if classification_type == ClassificationType.CDNA_REFERENCE_AGREE:
                    params["pos"] = int(params["pos"])
                    return CdnaReferenceAgreeClassification(**params)

                if classification_type == ClassificationType.CDNA_DELINS:
                    params["pos0"] = int(params["pos0"])
                    params["pos1"] = (
                        int(params["pos1"])
                        if params["pos1"] is not None
                        else params["pos1"]
                    )
                    return CdnaDelInsClassification(**params)

                if classification_type == ClassificationType.CDNA_DELETION:
                    params["pos0"] = int(params["pos0"])
                    params["pos1"] = (
                        int(params["pos1"])
                        if params["pos1"] is not None
                        else params["pos1"]
                    )
                    return CdnaDeletionClassification(**params)

                if classification_type == ClassificationType.CDNA_INSERTION:
                    params["pos0"] = int(params["pos0"])
                    params["pos1"] = (
                        int(params["pos1"])
                        if params["pos1"] is not None
                        else params["pos1"]
                    )
                    return CdnaInsertionClassification(**params)

        return None

    def _genomic_classification(
        self, token: HgvsToken, params: Dict
    ) -> Optional[Classification]:
        """Determine if hgvs token matches regex checks and return corresponding genomic
        classification if a match is found. Only checks against 'simple'
        duplication/deletions.

        :param token: hgvs token
        :param params: Base fields for a classification. This will get mutated if a
            match is found.
        :return: genomic classification if hgvs token matches regex checks. Else, `None`
        """
        for regex, _, classification_type in GENOMIC_REGEXPRS:
            match = self._regex_match(token.change, regex)

            if match:
                match_dict = match.groupdict()
                params.update(match_dict)

                if classification_type == ClassificationType.GENOMIC_SUBSTITUTION:
                    params["pos"] = int(params["pos"])
                    return GenomicSubstitutionClassification(**params)

                if classification_type == ClassificationType.GENOMIC_REFERENCE_AGREE:
                    params["pos"] = int(params["pos"])
                    return GenomicReferenceAgreeClassification(**params)

                if classification_type == ClassificationType.GENOMIC_DELINS:
                    params["pos0"] = int(params["pos0"])
                    params["pos1"] = (
                        int(params["pos1"])
                        if params["pos1"] is not None
                        else params["pos1"]
                    )
                    return GenomicDelInsClassification(**params)

                if classification_type == ClassificationType.GENOMIC_INSERTION:
                    params["pos0"] = int(params["pos0"])
                    params["pos1"] = (
                        int(params["pos1"])
                        if params["pos1"] is not None
                        else params["pos1"]
                    )
                    return GenomicInsertionClassification(**params)

                if classification_type == ClassificationType.GENOMIC_DELETION:
                    params["pos0"] = int(params["pos0"])
                    params["pos1"] = (
                        int(params["pos1"])
                        if params["pos1"] is not None
                        else params["pos1"]
                    )
                    return GenomicDeletionClassification(**params)

                if classification_type == ClassificationType.GENOMIC_DUPLICATION:
                    params["pos0"] = int(params["pos0"])
                    params["pos1"] = (
                        int(params["pos1"])
                        if params["pos1"] is not None
                        else params["pos1"]
                    )
                    return GenomicDuplicationClassification(**params)

        return None

    def _genomic_ambiguous_classification(
        self, token: HgvsToken, params: Dict
    ) -> Optional[Classification]:
        """Determine if hgvs token matches regex checks and return corresponding genomic
        ambiguous classification if a match is found. Only checks against ambiguous
        duplication/deletions.

        :param token: hgvs token
        :param params: Base fields for a classification. This will get mutated if a
            match is found.
        :return: genomic ambiguous classification if hgvs token matches regex checks.
            Else, `None`
        """
        if token.token.endswith("dup"):
            return self._genomic_dup_ambiguous_classification(token, params)

        if token.token.endswith("del"):
            return self._genomic_del_ambiguous_classification(token, params)

        return None

    @staticmethod
    def _update_ambiguous_params(params: Dict, regex_type: AmbiguousRegexType) -> None:
        """Mutates `params` to match correct types and gets associated ambiguous type
        from fields in `params`

        :param params: Fields for a classification. This will get mutated.
        :param regex_type: The kind of ambiguous regex that was used
        """
        params["pos0"] = (
            int(params["pos0"]) if params["pos0"] != "?" else params["pos0"]
        )

        if "pos1" in params:
            params["pos1"] = (
                int(params["pos1"]) if params["pos1"] != "?" else params["pos1"]
            )
        else:
            params["pos1"] = None

        params["pos2"] = (
            int(params["pos2"]) if params["pos2"] != "?" else params["pos2"]
        )

        if "pos3" in params:
            params["pos3"] = (
                int(params["pos3"]) if params["pos3"] != "?" else params["pos3"]
            )
        else:
            params["pos3"] = None

        ambiguous_type = get_ambiguous_type(
            params["pos0"], params["pos1"], params["pos2"], params["pos3"], regex_type
        )
        if ambiguous_type:
            params["ambiguous_type"] = ambiguous_type

    def _genomic_dup_ambiguous_classification(
        self, token: HgvsToken, params: Dict
    ) -> Optional[Classification]:
        """Determine if hgvs token matches regex checks and return corresponding genomic
        ambiguous duplication classification if a match is found. Only checks against
        genomic ambiguous duplications.

        :param token: hgvs token
        :param params: Base fields for a classification. This will get mutated if a
            match is found.
        :return: genomic ambiguous duplication classification if hgvs token matches
            regex checks. Else, `None`
        """
        for regex, _, classification_type, regex_type in GENOMIC_DUP_AMBIGUOUS_REGEXPRS:
            match = regex.match(token.change)

            if match:
                match_dict = match.groupdict()
                params.update(match_dict)

                if (
                    classification_type
                    == ClassificationType.GENOMIC_DUPLICATION_AMBIGUOUS
                ):
                    self._update_ambiguous_params(params, regex_type)

                    # If ambiguous type not in params, it means we don't support it yet
                    if "ambiguous_type" in params:
                        return GenomicDuplicationAmbiguousClassification(**params)
        return None

    def _genomic_del_ambiguous_classification(
        self, token: HgvsToken, params: Dict
    ) -> Optional[Classification]:
        """Determine if hgvs token matches regex checks and return corresponding genomic
        ambiguous deletion classification if a match is found. Only checks against
        genomic ambiguous deletion.

        :param token: hgvs token
        :param params: Base fields for a classification. This will get mutated if a
            match is found.
        :return: genomic ambiguous deletion classification if hgvs token matches regex
            checks. Else, `None`
        """
        for regex, _, classification_type, regex_type in GENOMIC_DEL_AMBIGUOUS_REGEXPRS:
            match = regex.match(token.change)

            if match:
                match_dict = match.groupdict()
                params.update(match_dict)

                if classification_type == ClassificationType.GENOMIC_DELETION_AMBIGUOUS:
                    self._update_ambiguous_params(params, regex_type)

                    # If ambiguous type not in params, it means we don't support it yet
                    if "ambiguous_type" in params:
                        return GenomicDeletionAmbiguousClassification(**params)
        return None
    


#
"""Module for to_vrs endpoint."""
from urllib.parse import unquote
from variation.classify import Classify
from variation.tokenize import Tokenize

class ToVRS:
    """The class for translating variation strings to VRS representations."""

    def __init__(
        self,
        #seqrepo_access: SeqRepoAccess,
        tokenizer: Tokenize,
        classifier: Classify,
        #validator: Validate,
        #translator: Translate,
    ) -> None:
        """Initialize the ToVRS class.

        :param SeqRepoAccess seqrepo_access: Access to SeqRepo
        :param Tokenize tokenizer: Tokenizer class for tokenizing
        :param Classify classifier: Classifier class for classifying tokens
        :param Validate validator: Validator class for validating valid inputs
        :param Translate translator: Translating valid inputs
        """
        #super().__init__(seqrepo_access)
        self.tokenizer = tokenizer
        self.classifier = classifier
        #self.validator = validator
        #self.translator = translator


   

    async def to_vrs(self, q: str) :
        """Return a VRS-like representation of all validated variations for a query.

        :param str q: The variation to translate (HGVS, gnomAD VCF, or free text) on
            GRCh37 or GRCh38 assembly
       # :return: ToVRSService containing VRS variations and warnings
        """
        warnings = []
        variations = []
        params = {
            "search_term": q,
            "warnings": warnings,
        }

        # Get tokens for input query
        tokens = self.tokenizer.perform(unquote(q.strip()), warnings)
        if warnings:
            params["warnings"] = warnings
            return ERRORMESSAGE(**params)

        # Get classification for list of tokens
        classification = self.classifier.perform(tokens)
        if not classification:
            params["warnings"] = [f"Unable to find classification for: {q}"]
            return ERRORMESSAGE(**params)


        params["warnings"] = warnings
        return ToVRSService(**params)  #actually return classification type to allow for error messages