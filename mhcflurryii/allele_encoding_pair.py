from .allele_encoding import AlleleEncoding


class AlleleEncodingPair(object):
    def __init__(
            self,
            alpha_allele_encoding,
            beta_allele_encoding):
        """
        """

        self.alpha_allele_encoding = alpha_allele_encoding
        self.beta_allele_encoding = beta_allele_encoding

    def from_pairs(self, allele_pairs):
        alpha_alleles = [a for (a, b) in allele_pairs]
        beta_alleles = [b for (a, b) in allele_pairs]
        return AlleleEncodingPair(
            AlleleEncoding(
                alpha_alleles,
                borrow_from=self.alpha_allele_encoding),
            AlleleEncoding(
                beta_alleles,
                borrow_from=self.beta_allele_encoding),
        )

    @property
    def allele_encodings(self):
        return [
            ("alpha", self.alpha_allele_encoding),
            ("beta", self.beta_allele_encoding)
        ]

    @property
    def allele_pairs(self):
        return [
            (a, b)
            for (a, b)
            in zip(
                self.alpha_allele_encoding.alleles,
                self.beta_allele_encoding.alleles)
        ]
