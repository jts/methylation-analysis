# Utility functions

# Split a string formated like key1=a;key2=b into a dictionary
def str2dict(s):
    return dict( (k, v) for k, v in (pair.split("=") for pair in s.split(";")) )

# These fucntions take in a string used in filenames/output and prettyifies them for figures and tables
def display_sample_name(s):
    fields = s.split("_")
    assert(fields[0] == "ecoli")
    return "\\emph{E. coli}" + " " + fields[1].upper()

def display_model(s):
    model_dict = {'t.006':'template', 'c.p1.006':'comp.pop1', 'c.p2.006':'comp.pop2'}
    return model_dict[s]

# From: http://stackoverflow.com/questions/3154460/python-human-readable-large-numbers
def display_number(n):
    suffix= ['','K','M','B' ]
    n = float(n)
    millidx = max(0,min(len(suffix)-1,
                    int(math.floor(0 if n == 0 else math.log10(abs(n))/3))))

    return '{:.0f}{}'.format(n / 10**(3 * millidx), suffix[millidx])

def display_treatment(s):
    return s.replace("_", "+").replace("MSssI", "M.SssI").replace("pcr", "PCR")

