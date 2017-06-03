COUNT_PATH=$1
COND_PATH=$2
SAVE_DIR=$3
ALL_GENE=$4

if [ $# -ne 3 ] && [ $# -ne 4 ]; then
  echo "Invalid Usage"
else
  echo "================================"
  echo "  differential expressed gene   "
  echo "================================"
  Rscript differentialExpressedGene.R $COUNT_PATH $COND_PATH $SAVE_DIR;

  echo "================================"
  echo "  count table to expression     "
  echo "================================"
  python countTable2exprs.py $COUNT_PATH $SAVE_DIR/raw_exprs.csv;
  python csv2tsv.py < $SAVE_DIR/raw_exprs.csv > $SAVE_DIR/raw_exprs.tsv;

  echo "================================"
  echo "        gene regulation         "
  echo "================================"
  python genies_script.py $SAVE_DIR/raw_exprs.tsv $SAVE_DIR;

  echo "================================"
  echo "            Main Part           "
  echo "================================"
  if [ $# -eq 3 ]; then
    python upstream_regulator.py $SAVE_DIR/raw_exprs.tsv $SAVE_DIR/VIM.csv $SAVE_DIR/observedRegulation.csv $SAVE_DIR
  elif [ $# -eq 4 ]; then
    python upstream_regulator.py $SAVE_DIR/raw_exprs.tsv $SAVE_DIR/VIM.csv $SAVE_DIR/observedRegulation.csv $SAVE_DIR $ALL_GENE
  fi

fi
