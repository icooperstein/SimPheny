from simpheny.scoring import combine_fixed_ebm, simpheny_score, confidence_tier

def test_udn_rpl13_first_match_score():
    p = combine_fixed_ebm(0.0001, 0.0002, "udn")
    score = simpheny_score(p)
    assert abs(score - 6.090576) < 0.01

def test_confidence():
    assert confidence_tier(4.5) == "High"
    assert confidence_tier(3.0) == "Medium"
    assert confidence_tier(2.0) == "Low"
