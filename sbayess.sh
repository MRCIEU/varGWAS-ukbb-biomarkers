PATH="$PATH":/mnt/storage/home/ml18692/apps/sbayess/gctb_2.03beta_Linux

# Calling SBayesR
gctb \
--sbayes R \
--ldm ./data/sbayess/ukbEURu_hm3_shrunk_sparse/ukbEURu_hm3_chr22_v3_50k.ldm.sparse \
--pi 0.95,0.02,0.02,0.01 \
--gamma 0.0,0.01,0.1,1 \
--gwas-summary data/alkaline_phosphatase.30610.0.0.ma \
--chain-length 10000 \
--burn-in 2000 \
--out-freq 10 \
--out chr22