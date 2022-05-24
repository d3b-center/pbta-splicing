cat  cc_bs_ids.txt | awk '{if($2~'6'){ print $1"\tC6"}}' > cc_bs_ids_c6.txt
cat  cc_bs_ids.txt | awk '{if($2~'5'){ print $1"\tC5"}}' > cc_bs_ids_c5.txt
cat  cc_bs_ids.txt | awk '{if($2~'4'){ print $1"\tC4"}}' > cc_bs_ids_c4.txt
cat  cc_bs_ids.txt | awk '{if($2~'3'){ print $1"\tC3"}}' > cc_bs_ids_c3.txt
cat  cc_bs_ids.txt | awk '{if($2~'2'){ print $1"\tC2"}}' > cc_bs_ids_c2.txt
cat  cc_bs_ids.txt | awk '{if($2~'1'){ print $1"\tC1"}}' > cc_bs_ids_c1.txt

cat cc_bs_ids_c1.txt cc_bs_ids_c2.txt cc_bs_ids_c3.txt cc_bs_ids_c4.txt cc_bs_ids_c5.txt | awk '{print $1"\tOther" }' > cc_bs_ids_notc6.txt
cat cc_bs_ids_c1.txt cc_bs_ids_c2.txt cc_bs_ids_c3.txt cc_bs_ids_c4.txt cc_bs_ids_c6.txt | awk '{print $1"\tOther" }' > cc_bs_ids_notc5.txt
cat cc_bs_ids_c1.txt cc_bs_ids_c2.txt cc_bs_ids_c3.txt cc_bs_ids_c5.txt cc_bs_ids_c6.txt | awk '{print $1"\tOther" }' > cc_bs_ids_notc4.txt
cat cc_bs_ids_c1.txt cc_bs_ids_c2.txt cc_bs_ids_c4.txt cc_bs_ids_c5.txt cc_bs_ids_c6.txt | awk '{print $1"\tOther" }' > cc_bs_ids_notc3.txt
cat cc_bs_ids_c1.txt cc_bs_ids_c3.txt cc_bs_ids_c4.txt cc_bs_ids_c5.txt cc_bs_ids_c6.txt | awk '{print $1"\tOther" }' > cc_bs_ids_notc2.txt
cat cc_bs_ids_c2.txt cc_bs_ids_c3.txt cc_bs_ids_c4.txt cc_bs_ids_c5.txt cc_bs_ids_c6.txt | awk '{print $1"\tOther" }' > cc_bs_ids_notc1.txt
 
