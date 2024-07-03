import typing as tp
from Bio.Seq import Seq
from Bio.SeqUtils import CheckSum

class SeqIDGeneratorException(Exception):
    pass

class SeqIDGenerator(object):
    """
    SeqIDGenerator generates unique identifiers for sequences using different encoding methods.
    
    Attributes:
        defined_encoders (list): List of supported encoders.
        default_encoder (function): The default encoder function to use if none is specified.
    """
    
    defined_encoders = [
        "crc32",
        "crc64",
        "gcg",
        "seguid"
    ]
    def __init__(
        self,
        default_encoder: tp.Optional[str]=None
    ):
        if default_encoder is not None:
            self.default_encoder = self._get_encoder_function(default_encoder)
        else:
            self.default_encoder = None
    def _get_encoder_function(
        self,
        encoder: str
    ):
        if not isinstance(encoder, str):
            raise TypeError(
                "'encoder' has to be a str"
            )
        if not encoder in self.defined_encoders:
            raise ValueError(
                "'encoder' accepts only one of following: "
                f"{self.defined_encoders}"
            )
        return getattr(CheckSum, encoder)
    def encode(
        self,
        seq: tp.Union[str, Seq],
        encoder: tp.Optional[str]=None
    ):
        if encoder is None:
            if not self.default_encoder:
                raise SeqIDGeneratorException(
                    "There is no default encoder defined; 'encoder' can not be None"
                )
            encoder = self.default_encoder
        else:
            encoder = self._get_encoder_function(encoder)
        return encoder(seq.upper())
