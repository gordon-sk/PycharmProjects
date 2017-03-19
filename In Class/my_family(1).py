__author__ = 'ghemingway'


my_family = {
    "graham":   {"mother": "anne", "father": "tim", "children": ["penny", "archer"], "year": 1976},
    "anne":     {"mother": "betty", "father": "david m", "children": ["graham", "kristen"], "year": 1948},
    "tim":      {"mother": "lucy", "father": "david h", "children": ["graham", "kristen"], "year": 1946},
    "kristen":  {"mother": "anne", "father": "tim", "children": ["dylan"], "year": 1973},
    "david h":  {"mother": "", "father": "", "children": ["tim", "stuart", "sam", "bee", "margaret"], "year": 1913},
    "david m":  {"mother": "", "father": "", "children": ["anne", "claire", "jane"], "year": 1907},
    "betty":    {"mother": "", "father": "", "children": ["anne", "claire", "jane"], "year": 1922},
    "claire":   {"mother": "betty", "father": "david m", "children": ["zabby", "willy"], "year": 1946},
    "jane":     {"mother": "betty", "father": "david m", "children": [], "year": 1951},
    "lucy":     {"mother": "", "father": "", "children": ["tim", "stuart", "sam", "bee", "margaret"], "year": 1914},
    "penny":    {"mother": "celeste", "father": "graham", "children": [], "year": 2009},
    "archer":   {"mother": "celeste", "father": "graham", "children": [], "year": 2013},
    "stuart":   {"mother": "lucy", "father": "david h", "children": [], "year": 1942},
    "sam":      {"mother": "lucy", "father": "david h", "children": [], "year": 1952},
    "billy":    {"mother": "", "father": "", "children": ["zabby", "willy"], "year": 1946},
    "bee":      {"mother": "lucy", "father": "david h", "children": [], "year": 1948},
    "margaret": {"mother": "lucy", "father": "david h", "children": [], "year": 1944},
    "celeste":  {"mother": "", "father": "", "children": ["penny", "archer"], "year": 1979},
    "dylan":    {"mother": "kristen", "father": "rob", "children": [], "year": 2006},
    "rob":      {"mother": "", "father": "", "children": ["dylan"], "year": 1969},
    "zabby":    {"mother": "claire", "father": "billy", "children": [], "year": 1983},
    "willy":    {"mother": "claire", "father": "billy", "children": [], "year": 1977},
}


def print_ancestors(family, person, t=0, year=1900):
    if family[person]['mother'] != '':
        print(family[person]['mother'] + ':\t' + str(t))
        print_ancestors(family, family[person]['mother'], t+1)
    if family[person]['father'] != '':
        print(family[person]['father'] + ':\t' + str(t))
        print_ancestors(family, family[person]['father'], t+1)


def print_descendants(family, person, t=0, year=2014):
    pass


if __name__ == "__main__":
    print_ancestors(my_family, "penny")
    print_descendants(my_family, "betty")
