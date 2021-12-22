# %%
import tkinter as tk
import pandas as pd

################################ Functions ################################

def calc_param(prot_seq, red_cyst):
    global mw
    global ext_coeff
    global spec_abs
    global aaa_count
    global naa
    global aliph_index
    global res_count
    global negative_aa
    global positive_aa
    global tot_neg_res
    global tot_pos_res
    global aa_list
    aa_mw = {"a":71.0779, "r":156.18568, "n":114.10264, "d":115.0874, "c":103.1429, "e":129.11398, "q":128.12922, "g":57.05132, 
            "h":137.13928, "i":113.15764, "l":113.15764, "k":128.17228, "m":131.19606, "f":147.17386, "p":97.11518, "s":87.0773, 
            "t":101.10388, "w":186.2099, "y":163.17326, "v":99.13106}
    letter_code = {"a":"Alanine", "r":"Arginine", "n":"Asparagine", "d":"Aspartic acid", "c":"Cysteine", "e":"Glutamic acid", "q":"Glutamine", "g":"Glutamine", 
                  "h":"Histidine", "i":"Isoleucine", "l":"Leucine", "k":"Lysine", "m":"Methionine", "f":"Phenylalanine", "p":"Proline", "s":"Serine", 
                  "t":"Threonine", "w":"Tryptophan", "y":"Tyrosine", "v":"Valine"}
    neg_aa = ["d", "e"]
    pos_aa = ["r", "k"]
    negative_aa = {}
    positive_aa = {}
    mw = 0
    naa = 0
    res_count = {}
    for char in prot_seq.lower():
            if char in aa_mw:
                mw += aa_mw[char]
                naa += 1
                if char in res_count:
                    res_count[char] += 1
                else:
                    res_count[char] = 1 
                if char in neg_aa:
                    if char in negative_aa:
                        negative_aa[char] += 1
                    else:
                        negative_aa[char] = 1
                if char in pos_aa:
                    if char in positive_aa:
                        positive_aa[char] += 1
                    else:
                        positive_aa[char] = 1

    mw += 18.0105
    
    aaa_count = {}
    aaa = ["y", "w", "c"]
    for char in prot_seq.lower():
        if char in aaa:
            if char in aaa_count:
                aaa_count[char] += 1
            else:
                aaa_count[char] = 1
        else:
            pass
    for char in aaa:
        if not char in aaa_count:
            aaa_count[char] = 0 
    if red_cyst:
        ext_coeff = aaa_count["y"]*1490 + aaa_count["w"]*5500
    else:
        ext_coeff = aaa_count["y"]*1490 + aaa_count["w"]*5500 + aaa_count["c"]*125
    spec_abs = ext_coeff / mw

    aliph_index = (res_count["a"]/naa + 2.9 * res_count["v"]/naa + 3.9 * (res_count["i"]/naa + res_count["l"]/naa)) * 100

    tot_neg_res = negative_aa["d"] + negative_aa["e"]
    tot_pos_res = positive_aa["r"] + positive_aa["k"]

    aa_list = {}
    for ch in prot_seq.lower():
        if ch in aa_mw:
            if ch in aa_list:
                pass
            else:
                aa_list[ch] = [letter_code[ch], res_count[ch], res_count[ch]/naa*100]
        else:
            pass

    return mw, ext_coeff, spec_abs, aaa_count, naa, aliph_index, tot_neg_res

def protein_parameters():
    prot_seq = protein_sequence.get("1.0", tk.END)
    red_cyst = boolvar.get()
    result_name = entry.get()
    calc_param(prot_seq, red_cyst)
    protein_composition = pd.DataFrame.from_dict(aa_list, orient='index', columns=["Amino Acid", "Number of Residues", "Fraction of Total Protein (%)"])
    protein_composition.to_excel('{}_protein_composition.xlsx'.format(result_name))

    params = """{} Protein Parameters

    1. Molecular weight (Da): {}
    2. Molar extinction coefficient (M-1 cm-1): {}
    3. Absorbance at 0.1 mg/mL: {}
    4. Protein length: {}
    5. Aromatic amino acids: {}
    6. Aliphatic index: {}
    7. Number of negatively charged residues: {}
    8. Number of positively charged residues: {}""".format(result_name.upper(), mw, ext_coeff, spec_abs, naa, aaa_count, aliph_index, tot_neg_res, tot_pos_res)
    f = open(result_name + ".txt", "x")
    f.write(params)
    f.close()

################################ Application ################################

# A Tk class object has to always be created at the beginning of the code
window = tk.Tk()

################################ First frame (logo) ################################

logo_frame = tk.Frame(master=window, borderwidth=1, padx=5, pady=5)
logo_frame.grid(row=0, column=0)

primary_font = ('Palatino Linotype', 10)
other_font = ("Palatino Linotype", 8)

logo = tk.Label(
        text="  Esteban's little protein parameter calculator  ",
        foreground='white',
        background='#505050',
        width=56,
        height=2,
        master=logo_frame,
        relief=tk.RIDGE,
        borderwidth=2
        )
# pack() tells tkinter to size the window as small as it can
logo.pack()

################################ Second frame (Protein name) ################################

protein_name_frame = tk.Frame(master=window, borderwidth=1)
protein_name_frame.grid(row=1, column=0, sticky="W", pady=5, padx=5)

protein_name = tk.Label(
        text='  Name of your protein  ',
        fg = 'white',
        background='#4682B4',
        # relief= tk.RAISED,
        width=18,
        height=2,
        master=protein_name_frame
)

protein_name.pack(side=tk.LEFT)

# test = tk.Text(master=window, background='white')
# test.pack()

entry = tk.Entry(
        fg='black',
        bg='white',
        width=35,
        master=protein_name_frame,
        borderwidth=1,
        font=other_font
        )

entry.pack(pady=5, padx=5)

################################ Third frame (Protein sequence) ################################

protein_sequence_frame = tk.Frame(master=window)
protein_sequence_frame.grid(row=2, column=0, sticky='nsew', padx=5)

protein_seq_label = tk.Label(
        text=' Insert your protein sequence below (one letter code) ',
        fg = 'white',
        background='#4682B4',
        # relief=tk.RAISED,
        borderwidth=2,
        width=41,
        height=2,
        master=protein_sequence_frame
)

protein_seq_label.grid(row=0, column=0, sticky='NW', pady=4)

protein_sequence = tk.Text(
        master=protein_sequence_frame,
        # relief=tk.RAISED,
        borderwidth=1,
        height=10,
        width=66,
        font=other_font
)
protein_sequence.grid(row=1, column=0, sticky='nsew')

################################ Fourth frame (reduced cysteines checkbutton) ################################

reduced_cysteines_frame = tk.Frame()

boolvar = tk.BooleanVar()
boolvar.set(False)

reduced_cysteines = tk.Checkbutton(
        master=reduced_cysteines_frame,
        text='All cysteines are reduced',
        variable=boolvar,
        onvalue=True
)

reduced_cysteines_frame.grid()

reduced_cysteines.pack()

################################ Fifth frame (calculate protein parameters button) ################################

button_frame = tk.Frame(relief=tk.RIDGE, pady=5)

button = tk.Button(
        text="Calculate protein parameters",
        width=25,
        height=1,
        bg='#4682B4',
        fg='white',
        master=button_frame,
        pady=5,
        relief=tk.RAISED,
        command=protein_parameters,
        borderwidth=2
        )

button.pack()
button_frame.grid()

# mainloop() runs the mainloop, which tell python to listen to events
window.mainloop()