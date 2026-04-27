import pytest

from prescriptions_processing import CodesExtractor

def test__tokenize_drug_name():
    assert CodesExtractor._tokenize_drug_name('Ala ma kota') == ' ala ma kota '
    assert CodesExtractor._tokenize_drug_name('!"#%&\'()*+,-./:;<=>?@[\\]^_`|~') is None
    assert CodesExtractor._tokenize_drug_name('1!"#%&\'()*+,-./:;<=>?@[\\]^_`|~a') == ' 1 a '
    assert CodesExtractor._tokenize_drug_name('new lines\nand tabs\t goes away') == ' new lines and tabs goes away '
    assert CodesExtractor._tokenize_drug_name('.7 dot.dot 9.0d 9.0 .9 9/.10 dot.9 9.dot .7') == ' 0.7 dot dot 9.0d 9.0 0.9 9 0.10 dot0.9 9 dot 0.7 '
    assert CodesExtractor._tokenize_drug_name('9.87.65 .8dot.8') == ' 9.87.65 0.8dot0.8 '
    assert CodesExtractor._tokenize_drug_name('1,234,567 12,3456 123,45ml Comma,comma comma,123') == ' 1234567 123456 123 45ml comma comma comma 123 '
    assert CodesExtractor._tokenize_drug_name('%1% %6 7   % 6 p%') == ' 1% 6 7% 6 p '
    assert CodesExtractor._tokenize_drug_name('HPV (Type 6,11,16,18)_Vac 0.5ml Pfs') == ' hpv type 6 11 16 18 vac 0.5ml pfs '
    assert CodesExtractor._tokenize_drug_name('Baxter_Chlorhex/Cetrimide .015%/.15% 1L') == ' baxter chlorhex cetrimide 0.015% 0.15% 1l '
    assert CodesExtractor._tokenize_drug_name('Hydrocortisone And Miconazole Cream 1 % + 2 %') == ' hydrocortisone and miconazole cream 1% 2% '

def test__remove_informative_suffix():
    assert CodesExtractor._remove_informative_suffix(' adren hcl ') == ' adren hcl '
    assert CodesExtractor._remove_informative_suffix(' adren crm ') == ' adren '
    assert CodesExtractor._remove_informative_suffix(" frumil systemic ") == ' frumil '
    assert CodesExtractor._remove_informative_suffix(" tryptizol cap inj tab ") == ' tryptizol '
    assert CodesExtractor._remove_informative_suffix(" valsartan tab amlodipine ") == " valsartan tab amlodipine "
    with pytest.raises(ValueError):
        CodesExtractor._remove_informative_suffix(" systemic tab ")
    with pytest.raises(TypeError):
        CodesExtractor._remove_informative_suffix(1234)
