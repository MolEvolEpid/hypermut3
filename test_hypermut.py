import pytest

from hypermut import *

def test_isfixedwidth():
  assert isfixedwidth('A') == True
  assert isfixedwidth('A|G') == True
  assert isfixedwidth('A|GG') == False

def test_die_widthnotfixed():
  with pytest.raises(SystemExit): 
    die_widthnotfixed('upstream') 

def test_widthonly():
  assert widthonly('A') == '.'
  assert widthonly('A|G') == '.'
  assert widthonly('AT') == '..'
  with pytest.raises(ValueError): 
    widthonly('A|GG')

def test_check_chars():
  check_chars('ABCD', list('ABCD'), 'Error')
  with pytest.raises(ValueError): 
    check_chars('ABCD', list('ABC'), 'Error')

def test_check_width():
  assert check_width('G,A,G,A,.,RD,.,YN|RC') == 2
  assert check_width('G,A,G,A,A,RD,B,YN|RC') == 3
  with pytest.raises(SystemExit): 
    check_width('G,A,G,A,.,RD,G|TT,YN|RC')
  with pytest.raises(SystemExit): 
    check_width('G,A,G,A,.,RD,.,Y|RC')
  with pytest.raises(SystemExit): 
    check_width('G|TT,A,G,A,.,RD,.,YN|RC')
  with pytest.raises(SystemExit): 
    check_width('G,A,G|TT,A,.,RD,.,YN|RC')
  with pytest.raises(ValueError): 
    check_width('G,A,G,A,.,R,.,YN|RC')
  with pytest.raises(ValueError): 
    check_width('G,A,G,A,A,RD,AG,YN|RC')

def test_check_context_overlap():
  check_context('G,A,G,A,.,RD,.,YN|RC')
  with pytest.raises(ValueError): 
    check_context('G,A,G,A,.,RD|RD,.,YN|RC') 
  with pytest.raises(ValueError): 
    check_context('G,A,G,A,.,RN,.,YN|RC')
  with pytest.warns(UserWarning):
    check_context('G,A,G,A,.,RD,.,YC')

def test_check_input_patterns():
  check_input_patterns('G,A,G,A,.,RD,.,YN|RC', list('GARDYNC'))
  with pytest.raises(ValueError): 
    check_input_patterns('Z,A,G,A,.,RD,.,YN|RC', list('GARDYNC'))

def test_allow_gaps_multistate():
  assert allow_gaps_multistate('A', 'strict') == '[-]*[A][-]*'
  assert allow_gaps_multistate('R', 'strict') == '[-]*[AGR][-]*'
  assert allow_gaps_multistate('AR', 'strict') == '[-]*[A][-]*[-]*[AGR][-]*'
  assert allow_gaps_multistate('A|R', 'strict') == '[-]*[A][-]*|[-]*[AGR][-]*'
  assert allow_gaps_multistate('A', 'partial') == '[-]*[ARDHVNWM][-]*'
  assert allow_gaps_multistate('R', 'partial') == '[-]*[AGRBDHVNWSKM][-]*'
  assert allow_gaps_multistate('AR', 'partial') == '[-]*[ARDHVNWM][-]*[-]*[AGRBDHVNWSKM][-]*'
  assert allow_gaps_multistate('A|R', 'partial') == '[-]*[ARDHVNWM][-]*|[-]*[AGRBDHVNWSKM][-]*'

def test_get_potential_matches():
  assert [(m.start(0), m.end(0)) for m in get_potential_matches('TGT', 0, None, re.compile('[-]*[G][-]*', re.IGNORECASE))] == [(1,2)]
  assert [(m.start(0), m.end(0)) for m in get_potential_matches('TRT', 0, None, re.compile('[-]*[AGR][-]*', re.IGNORECASE))] == [(1,2)]
  assert [(m.start(0), m.end(0)) for m in get_potential_matches('AGT', 0, None, re.compile('(?<=)[-]*([ACGTRYBDHVNWSKM])[-]*(?=[-]*[G][-]*[-]*[T][-]*)', re.IGNORECASE))] == [(0,1)]
  assert [(m.start(0), m.end(0)) for m in get_potential_matches('AGGG', 0, None, re.compile('(?<=)[-]*([ACGTRYBDHVNWSKM])[-]*(?=[-]*[G][-]*[-]*[G][-]*)', re.IGNORECASE))] == [(0, 1), (1, 2)]
  assert [(m.start(0), m.end(0)) for m in get_potential_matches('TGT', 0, 1, re.compile('[-]*[G][-]*', re.IGNORECASE))] == []
  assert [(m.start(0), m.end(0)) for m in get_potential_matches('AGT', 0, 2, re.compile('(?<=)[-]*([ACGTRYBDHVNWSKM])[-]*(?=[-]*[G][-]*[-]*[T][-]*)', re.IGNORECASE))] == []
  assert [(m.start(0), m.end(0)) for m in get_potential_matches('AGGG', 0, 3, re.compile('(?<=)[-]*([ACGTRYBDHVNWSKM])[-]*(?=[-]*[G][-]*[-]*[G][-]*)', re.IGNORECASE))] == [(0, 1)]

def test_compile_regex():
  assert compile_regex('G') == re.compile('G', re.IGNORECASE)
  assert compile_regex('G', 'A', 'T') == re.compile("(?<=A)[-]*(G)[-]*(?=T)",re.I)

iupac_dict = {"A": ["A"], "C": ["C"], "G": ["G"], "T": ["T"],
                "R": ["A", "G"], "Y": ["C", "T"], "S": ["G", "C"], "W": ["A", "T"], "K": ["G", "T"], "M": ["A", "C"],
                "B": ["C", "G", "T"], "D": ["A", "G", "T"], "H": ["A", "C", "T"], "V": ["A", "C", "G"],
                "N": ["A", "C", "G", "T"]}

def test_compute_context_prop():
  assert compute_context_prop('ACGT', ['.'], iupac_dict) == 1
  assert compute_context_prop('A', ['A'], iupac_dict) == 1
  assert compute_context_prop('A', ['R'], iupac_dict) == 1
  assert compute_context_prop('R', ['A'], iupac_dict) == 0.5
  assert compute_context_prop('R', ['A', 'G'], iupac_dict) == 1
  with pytest.raises(AssertionError): 
    assert compute_context_prop('ACGT', 'A', iupac_dict) == 1

def test_find_match_weight():
  assert find_match_weight('AGT', 0, 1, 'A', ['.'], ['RD'], iupac_dict, 'strict') == (1,1)
  assert find_match_weight('AGT', 0, 1, 'A', ['.'], ['RD'], iupac_dict, 'partial') == (1,1)
  assert find_match_weight('TGT', 0, 1, 'A', ['.'], ['RD'], iupac_dict, 'strict') == (1,0)
  assert find_match_weight('TGT', 0, 1, 'A', ['.'], ['RD'], iupac_dict, 'partial') == (1,0)
  assert find_match_weight('AGC', 0, 1, 'A', ['.'], ['RD'], iupac_dict, 'strict') == (0,0)
  assert find_match_weight('AGC', 0, 1, 'A', ['.'], ['RD'], iupac_dict, 'partial') == (0,0)
  assert find_match_weight('A-GC', 1, 2, 'A', ['.'], ['RD'], iupac_dict, 'strict') == (0,0)
  assert find_match_weight('A-GC', 1, 2, 'A', ['.'], ['RD'], iupac_dict, 'partial') == (0,0)
  assert find_match_weight('A-GT', 0, 1, 'A', ['.'], ['RD'], iupac_dict, 'strict') == (1,1)
  assert find_match_weight('A-GT', 0, 1, 'A', ['.'], ['RD'], iupac_dict, 'partial') == (1,1)
  assert find_match_weight('RGT', 0, 1, 'A', ['.'], ['RD'], iupac_dict, 'strict') == (0,0)
  assert find_match_weight('RGT', 0, 1, 'A', ['.'], ['RD'], iupac_dict, 'partial') == (1,0.5)
  assert find_match_weight('AGR', 0, 1, 'A', ['.'], ['RD'], iupac_dict, 'strict') == (1,1)
  assert find_match_weight('AGR', 0, 1, 'A', ['.'], ['RD'], iupac_dict, 'partial') == (1,1)
  assert find_match_weight('ANT', 0, 1, 'A', ['.'], ['RD'], iupac_dict, 'strict') == (0,0)
  assert find_match_weight('ANT', 0, 1, 'A', ['.'], ['RD'], iupac_dict, 'partial') == (0.5,0.5)
  assert find_match_weight('RNT', 0, 1, 'A', ['.'], ['RD'], iupac_dict, 'strict') == (0,0)
  assert find_match_weight('RNT', 0, 1, 'A', ['.'], ['RD'], iupac_dict, 'partial') == (0.5,0.25)
  assert find_match_weight('A-GN', 0, 1, 'A', ['.'], ['RD'], iupac_dict, 'strict') == (0,0)
  assert find_match_weight('A-GN', 0, 1, 'A', ['.'], ['RD'], iupac_dict, 'partial') == (0.75,0.75)
  assert find_match_weight('TAAGT', 2, 3, 'A', ['A'], ['RD'], iupac_dict, 'strict') == (1,1)
  assert find_match_weight('TAAGT', 2, 3, 'A', ['A'], ['RD'], iupac_dict, 'strict') == (1,1)
  assert find_match_weight('AAGT', 1, 2, 'A', ['A'], ['RD'], iupac_dict, 'strict') == (1,1)
  assert find_match_weight('ATGT', 1, 2, 'A', ['A'], ['RD'], iupac_dict, 'strict') == (1,0)
  assert find_match_weight('AAGC', 1, 2, 'A', ['A'], ['RD'], iupac_dict, 'strict') == (0,0)
  assert find_match_weight('A-GC', 1, 2, 'A', ['A'], ['RD'], iupac_dict, 'strict') == (0,0)
  assert find_match_weight('-AA-GT', 2, 3, 'A', ['A'], ['RD'], iupac_dict, 'strict') == (1,1)
  assert find_match_weight('ARGT', 1, 2, 'A', ['A'], ['RD'], iupac_dict, 'strict') == (0,0)
  assert find_match_weight('ARGT', 1, 2, 'A', ['A'], ['RD'], iupac_dict, 'partial') == (1,0.5)
  assert find_match_weight('AAGR', 1, 2, 'A', ['A'], ['RD'], iupac_dict, 'strict') == (1,1)
  assert find_match_weight('AANT', 1, 2, 'A', ['A'], ['RD'], iupac_dict, 'strict') == (0,0)
  assert find_match_weight('AANT', 1, 2, 'A', ['A'], ['RD'], iupac_dict, 'partial') == (0.5,0.5)
  assert find_match_weight('ARNT', 1, 2, 'A', ['A'], ['RD'], iupac_dict, 'strict') == (0,0)
  assert find_match_weight('ARNT', 1, 2, 'A', ['A'], ['RD'], iupac_dict, 'partial') == (0.5,0.25)
  assert find_match_weight('RAAT', 1, 2, 'A', ['A'], ['RD'], iupac_dict, 'partial') == (0.5,0.5)
  assert find_match_weight('RANT', 1, 2, 'A', ['A'], ['RD'], iupac_dict, 'partial') == (0.25,0.25)
  assert find_match_weight('RRNT', 1, 2, 'A', ['A'], ['RD'], iupac_dict, 'partial') == (0.25,0.5*0.25)
  assert find_match_weight('R-A-GD', 2, 3, 'A', ['A'], ['RD'], iupac_dict, 'partial') == (0.5,0.5)
  assert find_match_weight('AAGT', 1, 2, 'A', ['A'], ['RD', 'AC'], iupac_dict, 'strict') == (1,1)
  assert find_match_weight('ATGT', 1, 2, 'A', ['A'], ['RD', 'AC'], iupac_dict, 'strict') == (1,0)
  assert find_match_weight('AAAC', 1, 2, 'A', ['A'], ['RD', 'AC'], iupac_dict, 'strict') == (1,1)

def test_summarize_matches():
  assert summarize_matches('GGT', 'AGT', 0, None,
                           re.compile('G'), re.compile("(?<=)[-]*(A)[-]*(?=GT)",re.I),
                           'A', ['.'], ['RD'],
                           iupac_dict, 'strict', 1, 'test') == (1,1)
  
def test_calc_fisher():
  assert calc_fisher(69,4,52,1)[1] == 0.28266930809303686

