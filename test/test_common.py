from mhc2flurry.common import make_allele_pairs


def test_allele_pairs():
    alleles = [
        "HLA-DRB1*07:01",
        "HLA-DRB1*16:01",
        "HLA-DRB4*01:03",
        "HLA-DRB5*02:02",
        "HLA-DPA1*01:03",
        "HLA-DPB1*02:01",
        "HLA-DPB1*23:01",
        "HLA-DQA1*01:02",
        "HLA-DQA1*02:01",
        "HLA-DQB1*02:02",
        "HLA-DQB1*05:02",
    ]
    result = make_allele_pairs(alleles)

    assert result == [
        'HLA-DRA*01:01-DRB1*07:01',
        'HLA-DRA*01:01-DRB1*16:01',
        'HLA-DRA*01:01-DRB4*01:03',
        'HLA-DRA*01:01-DRB5*02:02',
        'HLA-DPA1*01:03-DPB1*02:01',
        'HLA-DPA1*01:03-DPB1*23:01',
        'HLA-DQA1*01:02-DQB1*02:02',
        'HLA-DQA1*01:02-DQB1*05:02',
        'HLA-DQA1*02:01-DQB1*02:02',
        'HLA-DQA1*02:01-DQB1*05:02',
    ]
