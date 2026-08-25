import pandas as pd
from simpheny.ranking import add_gene_rankings

def test_rpl13_style_top_two_aggregation():
    df = pd.DataFrame({
        "Reference_Source":["udn"]*5,
        "Reference_ID":["A","B","C","D","E"],
        "Gene_Hit":["RPL13","RPL13","RPL13","RYR1","RYR1"],
        "SimPheny_Score":[6.090576,5.245944,4.452706,2.165453,1.549537],
    })
    _, genes = add_gene_rankings(df)
    rpl13 = genes[genes["Gene_Hit"]=="RPL13"].iloc[0]
    ryr1 = genes[genes["Gene_Hit"]=="RYR1"].iloc[0]
    assert abs(rpl13["SimPheny_Gene_Score"] - 5.66826) < 1e-5
    assert abs(ryr1["SimPheny_Gene_Score"] - 1.857495) < 1e-5
    assert int(rpl13["Candidate_Rank"]) == 1
