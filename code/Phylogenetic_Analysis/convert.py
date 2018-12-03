import csv


def is_valid_loci(loci):
    if isinstance(loci, str) and len(loci) == 2 and loci[0] in "AGTC" and loci[1] in "AGTC":
        return True
    else:
        return False


ori_csv_file = open("sample.csv", "r")
csv_reader = csv.reader(ori_csv_file, delimiter=',')

rst_cvs_file = open("result.csv", "w", newline='')
csv_writer = csv.writer(rst_cvs_file, delimiter=',')

replace_dict_A = {"AA": "0", "AG": "1", "AT": "1", "AC": "1",
                  "GA": "2", "GG": "2", "GT": "2", "GC": "2",
                  "TA": "2", "TG": "2", "TT": "2", "TC": "2",
                  "CA": "2", "CG": "2", "CT": "2", "CC": "2", "": "?"}
replace_dict_G = {"AA": "2", "AG": "2", "AT": "2", "AC": "2",
                  "GA": "1", "GG": "0", "GT": "1", "GC": "1",
                  "TA": "2", "TG": "2", "TT": "2", "TC": "2",
                  "CA": "2", "CG": "2", "CT": "2", "CC": "2", "": "?"}
replace_dict_T = {"AA": "2", "AG": "2", "AT": "2", "AC": "2",
                  "GA": "2", "GG": "2", "GT": "2", "GC": "2",
                  "TA": "1", "TG": "1", "TT": "0", "TC": "1",
                  "CA": "2", "CG": "2", "CT": "2", "CC": "2", "": "?"}
replace_dict_C = {"AA": "2", "AG": "2", "AT": "2", "AC": "2",
                  "GA": "2", "GG": "2", "GT": "2", "GC": "2",
                  "TA": "2", "TG": "2", "TT": "2", "TC": "2",
                  "CA": "1", "CG": "1", "CT": "1", "CC": "0", "": "?"}

dict_choice_map = {"A": replace_dict_A, "G": replace_dict_G, "T": replace_dict_T, "C": replace_dict_C}


dict_choice = []
fst_line = True
lines = []
for row in csv_reader:
    lines.append(row)
    if fst_line:
        fst_line = False
        dict_choice = [""] * (len(row)-2)
    else:
        row = row[2:]
        for col_cnt in range(len(row)):
            loci = row[col_cnt]
            if is_valid_loci(loci):
                dict_choice[col_cnt] = loci[0]

fst_line = True
for row in lines:
    if fst_line:
        fst_line = False
        csv_writer.writerow(row)
    else:
        rst_row = row[0:2]
        row = row[2:]
        for col_cnt in range(len(row)):
            loci = row[col_cnt]
            if dict_choice[col_cnt] in dict_choice_map:
                cur_line_dict = dict_choice_map[dict_choice[col_cnt]]
                rst_row += cur_line_dict[loci]
            else:
                rst_row.append("")
        csv_writer.writerow(rst_row)
