python countTable2exprs.py raw_countTable.csv raw_exprs.csv
python csv2tsv.py < raw_exprs.csv > raw_exprs.tsv
cp raw_exprs.tsv unknown_expression_data.tsv
cp unknown_expression_data.tsv TIGRESS-PKG-2.1/data/
