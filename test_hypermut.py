import pytest

from hypermut import *

def test_check_chars():
  check_chars('ABCD', list('ABCD'), 'Error')
  with pytest.raises(ValueError): 
    check_chars('ABCD', list('ABC'), 'Error')

def test_check_width():
  check_width('RD', 'test')
  check_width('A', 'test')
  with pytest.raises(ValueError): 
    check_width('G|TT', 'test')
  with pytest.raises(ValueError): 
    check_width('GT|T', 'test')

def test_check_input_patterns():
  check_input_patterns('G','A','.','RD', iupac_dict)
  check_input_patterns('G','A','A','RD', iupac_dict)
  check_input_patterns('G','A','AR','RD', iupac_dict)
  with pytest.raises(ValueError): 
    check_input_patterns('Z','A','.','RD', iupac_dict)
  with pytest.raises(ValueError): 
    check_input_patterns('G','A','.','C|RD', iupac_dict)
  with pytest.raises(ValueError): 
    check_input_patterns('GT','A','.','RD', iupac_dict)
  with pytest.raises(ValueError): 
    check_input_patterns('G','A','.','RD|RD', iupac_dict)

# def test_get_potential_matches():
#   assert [(m.start(0), m.end(0)) for m in get_potential_matches('TGT', 0, None, re.compile('[-]*[G][-]*', re.IGNORECASE))] == [(1,2)]
#   assert [(m.start(0), m.end(0)) for m in get_potential_matches('TRT', 0, None, re.compile('[-]*[AGR][-]*', re.IGNORECASE))] == [(1,2)]
#   assert [(m.start(0), m.end(0)) for m in get_potential_matches('AGT', 0, None, re.compile('(?<=)[-]*([ACGTRYBDHVNWSKM])[-]*(?=[-]*[G][-]*[-]*[T][-]*)', re.IGNORECASE))] == [(0,1)]
#   assert [(m.start(0), m.end(0)) for m in get_potential_matches('AGGG', 0, None, re.compile('(?<=)[-]*([ACGTRYBDHVNWSKM])[-]*(?=[-]*[G][-]*[-]*[G][-]*)', re.IGNORECASE))] == [(0, 1), (1, 2)]
#   assert [(m.start(0), m.end(0)) for m in get_potential_matches('TGT', 0, 1, re.compile('[-]*[G][-]*', re.IGNORECASE))] == []
#   assert [(m.start(0), m.end(0)) for m in get_potential_matches('AGT', 0, 2, re.compile('(?<=)[-]*([ACGTRYBDHVNWSKM])[-]*(?=[-]*[G][-]*[-]*[T][-]*)', re.IGNORECASE))] == []
#   assert [(m.start(0), m.end(0)) for m in get_potential_matches('AGGG', 0, 3, re.compile('(?<=)[-]*([ACGTRYBDHVNWSKM])[-]*(?=[-]*[G][-]*[-]*[G][-]*)', re.IGNORECASE))] == [(0, 1)]

iupac_dict = {"A": ["A"], "C": ["C"], "G": ["G"], "T": ["T"],
                "R": ["A", "G"], "Y": ["C", "T"], "S": ["G", "C"], "W": ["A", "T"], "K": ["G", "T"], "M": ["A", "C"],
                "B": ["C", "G", "T"], "D": ["A", "G", "T"], "H": ["A", "C", "T"], "V": ["A", "C", "G"],
                "N": ["A", "C", "G", "T"], "-": ["-"]}

def test_compute_context_prop():
  assert compute_context_prop('GCGT', 'ACGT', ['.'], 'D', iupac_dict) == 1
  assert compute_context_prop('GCGT', 'ACGT', ['.'], 'A', iupac_dict) == 1
  assert compute_context_prop('GCGT', 'ACGT', ['.'], 'B', iupac_dict) == 1
  assert compute_context_prop('G', 'A', ['A'], 'D', iupac_dict) == 1
  assert compute_context_prop('G', 'A', ['A'], 'A', iupac_dict) == 0
  assert compute_context_prop('A', 'A', ['A'], 'A', iupac_dict) == 1
  assert compute_context_prop('G', 'A', ['A'], 'B', iupac_dict) == 0.5
  assert compute_context_prop('A', 'A', ['A'], 'B', iupac_dict) == 1
  assert compute_context_prop('G', 'A', ['R'], 'D', iupac_dict) == 1
  assert compute_context_prop('G', 'R', ['A'], 'D', iupac_dict) == 0.5
  assert compute_context_prop('G', 'R', ['A', 'G'], 'D', iupac_dict) == 1
  assert compute_context_prop('GG', 'NT', ['RD'], 'D', iupac_dict) == 0.5
  assert compute_context_prop('GG', 'NT', ['YN', 'RC'], 'D', iupac_dict) == 0.5
  assert compute_context_prop('GG', 'WS', ['RD'], 'D', iupac_dict) == 0.25
  assert compute_context_prop('GG', 'WS', ['YN', 'RC'], 'D', iupac_dict) == 0.75
  assert compute_context_prop('GG', 'N', ['R'], 'D', iupac_dict) == 0.5
  assert compute_context_prop('GG', 'N', ['Y'], 'D', iupac_dict) == 0.5
  assert compute_context_prop('GG', 'NN', ['RR'], 'D', iupac_dict) == 0.25
  assert compute_context_prop('GG', 'NN', ['YY', 'RY', 'YR'], 'D', iupac_dict) == 0.75

def test_find_match_weight():
  assert find_match_weight('GGG', 'AGT', 0, 1, 'A', ['.'], ['RD'], 'D', iupac_dict, 'strict', False) == (1,1,0,0)
  assert find_match_weight('GGG', 'AGT', 0, 1, 'A', ['.'], ['RD'], 'D', iupac_dict, 'strict', True) == (1,1,0,0)
  assert find_match_weight('GGT', 'AGT', 0, 1, 'A', ['.'], ['RD'], 'A', iupac_dict, 'strict', False) == (1,1,0,0)
  assert find_match_weight('GGT', 'AGT', 0, 1, 'A', ['.'], ['RD'], 'B', iupac_dict, 'strict', False) == (1,1,0,0)
  assert find_match_weight('GGG', 'AGT', 0, 1, 'A', ['.'], ['RD'], 'D', iupac_dict, 'partial', False) == (1,1,0,0)
  assert find_match_weight('GGG', 'TGT', 0, 1, 'A', ['.'], ['RD'], 'D', iupac_dict, 'strict', False) == (1,0,0,0)
  assert find_match_weight('GGG', 'TGT', 0, 1, 'A', ['.'], ['RD'], 'D', iupac_dict, 'partial', False) == (1,0,0,0)
  assert find_match_weight('GGG', 'AGC', 0, 1, 'A', ['.'], ['RD'], 'D', iupac_dict, 'strict', False) == (0,0,1,1)
  assert find_match_weight('GGG', 'AGC', 0, 1, 'A', ['.'], ['RD'], 'D', iupac_dict, 'partial', False) == (0,0,1,1)
  assert find_match_weight('GGGG', 'A-GC', 1, 2, 'A', ['.'], ['RD'], 'D', iupac_dict, 'strict', False) == (0,0,0,0)
  assert find_match_weight('GGGG', 'A-GC', 1, 2, 'A', ['.'], ['RD'], 'D', iupac_dict, 'strict', True) == (0,0,0,0)
  assert find_match_weight('GGGG', 'A-GC', 1, 2, 'A', ['.'], ['RD'], 'D', iupac_dict, 'partial', False) == (0,0,0,0)
  assert find_match_weight('GGGG', 'A-GT', 0, 1, 'A', ['.'], ['RD'], 'D', iupac_dict, 'strict', False) == (1,1,0,0)
  assert find_match_weight('GGGG', 'A-GT', 0, 1, 'A', ['.'], ['RD'], 'D', iupac_dict, 'strict', True) == (0,0,0,0)
  assert find_match_weight('GGGG', 'A-GT', 0, 1, 'A', ['.'], ['RD'], 'D', iupac_dict, 'partial', False) == (1,1,0,0)
  assert find_match_weight('GGG', 'RGT', 0, 1, 'A', ['.'], ['RD'], 'D', iupac_dict, 'strict', False) == (0,0,0,0)
  assert find_match_weight('GGG', 'RGT', 0, 1, 'A', ['.'], ['RD'], 'D', iupac_dict, 'partial', False) == (1,0.5,0,0)
  assert find_match_weight('GGT', 'RGT', 0, 1, 'A', ['.'], ['RD'], 'A', iupac_dict, 'partial', False) == (1,0.5,0,0)
  assert find_match_weight('GGT', 'RGT', 0, 1, 'A', ['.'], ['RD'], 'B', iupac_dict, 'partial', False) == (1,0.5,0,0)
  assert find_match_weight('GGG', 'AGR', 0, 1, 'A', ['.'], ['RD'], 'D', iupac_dict, 'strict', False) == (1,1,0,0)
  assert find_match_weight('GGG', 'AGR', 0, 1, 'A', ['.'], ['RD'], 'A', iupac_dict, 'strict', False) == (1,1,0,0)
  assert find_match_weight('GGG', 'AGC', 0, 1, 'A', ['.'], ['RD'], 'B', iupac_dict, 'strict', False) == (0,0,0,0)
  assert find_match_weight('GGG', 'AGR', 0, 1, 'A', ['.'], ['RD'], 'D', iupac_dict, 'partial', False) == (1,1,0,0)
  assert find_match_weight('GGG', 'AGR', 0, 1, 'A', ['.'], ['RD'], 'A', iupac_dict, 'partial', False) == (1,1,0,0)
  assert find_match_weight('GGG', 'AGC', 0, 1, 'A', ['.'], ['RD'], 'B', iupac_dict, 'partial', False) == (0.5,0.5,0.5,0.5) # IS THIS WHAT WE WANT?
  assert find_match_weight('GGG', 'ANT', 0, 1, 'A', ['.'], ['RD'], 'D', iupac_dict, 'strict', False) == (0,0,0,0)
  assert find_match_weight('GGG', 'ANT', 0, 1, 'A', ['.'], ['RD'], 'D', iupac_dict, 'partial', False) == (0.5,0.5,0.5,0.5)
  assert find_match_weight('GGG', 'RNT', 0, 1, 'A', ['.'], ['RD'], 'D', iupac_dict, 'strict', False) == (0,0,0,0)
  assert find_match_weight('GGG', 'RNT', 0, 1, 'A', ['.'], ['RD'], 'D', iupac_dict, 'partial', False) == (0.5,0.25,0.5,0.25)
  assert find_match_weight('GGGG', 'A-GN', 0, 1, 'A', ['.'], ['RD'], 'D', iupac_dict, 'strict', False) == (0,0,0,0)
  assert find_match_weight('GGGG', 'A-GN', 0, 1, 'A', ['.'], ['RD'], 'D', iupac_dict, 'partial', False) == (0.75,0.75,0.25,0.25)
  assert find_match_weight('GGGGG', 'TAAGT', 2, 3, 'A', ['A'], ['RD'], 'D', iupac_dict, 'strict', False) == (1,1,0,0)
  assert find_match_weight('GGGGG', 'TAAGT', 2, 3, 'A', ['A'], ['RD'], 'D', iupac_dict, 'strict', False) == (1,1,0,0)
  assert find_match_weight('GGGG', 'AAGT', 1, 2, 'A', ['A'], ['RD'], 'D', iupac_dict, 'strict', False) == (1,1,0,0)
  assert find_match_weight('GGGG', 'ATGT', 1, 2, 'A', ['A'], ['RD'], 'D', iupac_dict, 'strict', False) == (1,0,0,0)
  assert find_match_weight('GGGG', 'AAGC', 1, 2, 'A', ['A'], ['RD'], 'D', iupac_dict, 'strict', False) == (0,0,1,1)
  assert find_match_weight('GGGG', 'A-GC', 1, 2, 'A', ['A'], ['RD'], 'D', iupac_dict, 'strict', False) == (0,0,0,0)
  assert find_match_weight('GGGGGG', '-AA-GT', 2, 3, 'A', ['A'], ['RD'], 'D', iupac_dict, 'strict', False) == (1,1,0,0)
  assert find_match_weight('GGGG', 'ARGT', 1, 2, 'A', ['A'], ['RD'], 'D', iupac_dict, 'strict', False) == (0,0,0,0)
  assert find_match_weight('GGGG', 'ARGT', 1, 2, 'A', ['A'], ['RD'], 'D', iupac_dict, 'partial', False) == (1,0.5,0,0)
  assert find_match_weight('GGGG', 'AAGR', 1, 2, 'A', ['A'], ['RD'], 'D', iupac_dict, 'strict', False) == (1,1,0,0)
  assert find_match_weight('GGGG', 'AANT', 1, 2, 'A', ['A'], ['RD'], 'D', iupac_dict, 'strict', False) == (0,0,0,0)
  assert find_match_weight('GGGG', 'AANT', 1, 2, 'A', ['A'], ['RD'], 'D', iupac_dict, 'partial', False) == (0.5,0.5,0.5,0.5)
  assert find_match_weight('GGGG', 'ARNT', 1, 2, 'A', ['A'], ['RD'], 'D', iupac_dict, 'strict', False) == (0,0,0,0)
  assert find_match_weight('GGGG', 'ARNT', 1, 2, 'A', ['A'], ['RD'], 'D', iupac_dict, 'partial', False) == (0.5,0.25,0.5,0.25)
  assert find_match_weight('GGGG', 'RAAT', 1, 2, 'A', ['A'], ['RD'], 'D', iupac_dict, 'partial', False) == (0.5,0.5,0.5,0.5)
  assert find_match_weight('GGGG', 'RANT', 1, 2, 'A', ['A'], ['RD'], 'D', iupac_dict, 'partial', False) == (0.25,0.25,0.75,0.75)
  assert find_match_weight('GGGG', 'RRNT', 1, 2, 'A', ['A'], ['RD'], 'D', iupac_dict, 'partial', False) == (0.25,0.5*0.25,0.75,0.5*0.75)
  assert find_match_weight('GGGGGG', 'R-A-GD', 2, 3, 'A', ['A'], ['RD'], 'D', iupac_dict, 'partial', False) == (0.5,0.5,0.5,0.5)
  assert find_match_weight('GGGG', 'AAGT', 1, 2, 'A', ['A'], ['RD', 'AC'], 'D', iupac_dict, 'strict', False) == (1,1,0,0)
  #assert find_match_weight('AAGT', 1, 2, 'A', ['A'], ['RD', 'RD'], iupac_dict, 'strict', False) == (1,1,0,0) # wrong but incorrect context caught earlier in script
  assert find_match_weight('GGGG', 'ATGT', 1, 2, 'A', ['A'], ['RD', 'AC'], 'D', iupac_dict, 'strict', False) == (1,0,0,0)
  assert find_match_weight('GGGG', 'AAAC', 1, 2, 'A', ['A'], ['RD', 'AC'], 'D', iupac_dict, 'strict', False) == (1,1,0,0)

  assert find_match_weight('GGG', 'ANN', 0, 1, 'A', ['.'], ['RR'], 'D', iupac_dict, 'partial', False) == (0.25,0.25,0.75,0.75)
  assert find_match_weight('GGG', 'ANN', 0, 1, 'A', ['.'], ['YY', 'YR', 'RY'], 'D', iupac_dict, 'partial', False) == (0.75,0.75,0.25,0.25)
  assert find_match_weight('GGG', 'NAN', 1, 2, 'A', ['R'], ['R'], 'D', iupac_dict, 'partial', False) == (0.25,0.25,0.75,0.75)
  assert find_match_weight('GGG', 'NAN', 1, 2, 'A', ['Y'], ['Y'], 'D', iupac_dict, 'partial', False) == (0.25,0.25,0.75,0.75)

def test_summarize_matches():
  assert summarize_matches('GGTT', 'AGTT', 0, None,
                           re.compile('G'), 'D',
                           'A', ['.'], ['RD'],
                           iupac_dict, 'strict', False, 1, 'test', None) == (1,1,1,0)
  assert summarize_matches('GAGAT', 'A-G-T', 0, None,
                           re.compile('G'), 'D',
                           'A', ['.'], ['RD'],
                           iupac_dict, 'strict', False, 1, 'test', None) == (1,1,0,0)
  
def test_calc_fisher():
  assert calc_fisher(69,4,52,1)[1] == 0.28266930809303686

