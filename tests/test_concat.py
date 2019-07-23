import os
from pandas_filter import concat




# ###################################################################################################

workdir_its = "data/comp_data/test_its"
workdir_ets = "data/comp_data/test_ets"

workdir_its = "data/Senecio/its"
workdir_ets = "data/Senecio/ets"

email = "martha.kandziora@yahoo.com"
percentage = 0.4
workdir_comb = "./tests/output/concat_test"

genelist = {"ITS": workdir_its,
            "ETS": workdir_ets}


conc = concat.Concat(workdir_comb, email)
#assert isinstance(conc, Concat)
conc = conc.run(genelistdict=genelist, workdir_comb=workdir_comb,
                       email=email, percentage=percentage, user_concat_fn=None)




