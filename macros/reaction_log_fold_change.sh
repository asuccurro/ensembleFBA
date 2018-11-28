
source ../venvpy/bin/activate
echo "" > ../outputs/numbers_reaction_log_fold_change.md
for r in Root9 Root491 Root66D1; do
    python pathwayAnalysis.py -o $r -e '_exclude_5N_and_nf_and_G12' -p '../outputs/' -f 'ensemble_50_size_26_gcs_11_ngcs_stochasticWeights_1' >> ../outputs/numbers_reaction_log_fold_change.md
    echo "">> ../outputs/numbers_reaction_log_fold_change.md
done

deactivate
