import doctest
import scipy.special
import math
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import os.path
import importlib

#####################
# Custom Exceptions #
#####################
class LambdaError(Exception):
    pass

class LambdaSyntaxError(LambdaError):
    pass

class LambdaEvaluationError(LambdaError):
    pass

######################################################
# Helper functions: Arithmetic in the lambda algebra #
######################################################
def tokenize(string):
    pass

def binom(n,m):
    '''
    Calculates binomial coefficients n choose m, with the convention that if n choose m = 0
    if n < m or m < 0
    '''
    # if n < m or m < 0:
    #     return 0
    return scipy.special.comb(n, m, exact=True)

def isadmiss(lst):
    '''
    helper function; checks to see if a tuple of numbers representing a monomial in the
    lambda algebra is admissible. 
    Outputs (True, None) if input is admissible, and outputs (False, i) if input is not admissible, 
    where i is the first index where admissibility fails.
    '''
    if not isinstance(lst, (list, tuple)):
        raise TypeError('Tried to check if nonlist is admissible') 
    for i in range(len(lst)-1):
        if 2*lst[i] >= lst[i + 1]:
            continue
        else:
            return (False, i)
    return (True, None)

def lexicographicorder(left, right):
    '''
    helper function; takes in two tuples of numbers, representing monomials in the lambda algebra,
    as inputs and checks to see if the left is less than the right. Return True if it is, False otherwise.
    '''
    if not isinstance(left, (list, tuple)) or not isinstance(right, (list, tuple)):
        raise TypeError("Tried to compare nonlists.")
    if len(left) == 0:
        if len(right) == 0:
            return False
        return True
    if len(right) == 0:
        return False

    if left[0] < right[0]:
        return True
    if left[0] > right[0]:
        return False
    return lexicographicorder(left[1:], right[1:])

def lexicographicorder2(left, right):
    '''
    helper function; takes in two lists of tuples of numbers, representing polynomials in the lambda algebra, as inputs
    and checks to see if the left is less than the right. Returns True is it is, False otherwise.
    '''
    if len(left) == 0:
        if len(right) == 0:
            return False
        return True
        
    if len(left) == 0:
        return False

    if lexicographicorder(left[0], right[0]):
        return True
    
    if left[0] == right[0]:
        return lexicographicorder2(left[1:], right[1:])

    return False

def mul(self, other):
    '''
    Helper function; returns 
    '''
    sign = []
    for mon1 in self:
        for mon2 in other:
            sign.append(mon1 + mon2)
    return sign

def mergesort(list, order):
    '''
    generic mergesort of a list, according to a specified order.
    '''
    if len(list) <= 1:
        return list
    left = mergesort(list[:len(list) // 2], order)
    right = mergesort(list[len(list) // 2:], order)

    prod = []
    leftplace = 0
    rightplace = 0

    while len(prod) != len(list):
        try:
            if order(left[leftplace], right[rightplace]):
                prod.append(left[leftplace])
                leftplace += 1
                continue
            prod.append(right[rightplace])
            rightplace += 1
        except:
            if len(left) == leftplace:
                prod = prod + right[rightplace:]
            else:
                prod = prod + left[leftplace:]

    return prod

def reduce(poly):
    dic = {}

    for mon in poly:
        if mon in dic:
            if dic[mon] == 1:
                dic[mon] = 0
            else:
                dic[mon] = 1
            continue
        dic[mon] = 1

    newpoly = []

    for i in dic.keys():
        if dic[i] % 2 == 1:
            newpoly.append(i)

    return mergesort(newpoly, lambda left, right : not lexicographicorder(left, right) and not left == right)

def admiss(poly):
    '''
    Takes in tuple representing monomial and makes it admissible. Returns list representing admissed polynomial.
    '''
    if len(poly) == 1:
        isit, ind = isadmiss(poly[0])
        if isit:
            return poly
        
        toosmall = poly[0][ind]
        toobig = poly[0][ind + 1]

        k = toobig - 2*toosmall
        
        if k <= 0:
            raise Exception('something weird happened')
        
        prod = []

        for j in range(k): #this might have to be k - 1; check later.
            if binom(k - 1 - j - 1, j) % 2 == 0:
                continue
            #good chance of arithmetic errors, check against Tangora later.
            prod += [(toosmall + k - 1 - j, 2*toosmall + 1 + j)]
        
        prod = mul(mul([poly[0][:ind]], prod), [poly[0][ind + 2:]])
        
        return admiss(prod)

    prod = []

    for mon in poly:
        prod += admiss([mon])
    
    return reduce(prod)

def bidegree(monomial):
    '''
    Returns bidegree of a monomial, in the form of a tuple (t - s, s). The bidegree of a generator \lambda_r is (t - s, s) = (r, 1)

    Inputs:
    monomial: a monomial in the lambda algebra, represented as a list of numbers
    '''
    return sum(monomial), len(monomial)

def ishomog(signature):
    '''
    checks a list/tuple of list/tuples to see if it is homogenous.
    '''
    if len(signature) <= 1:
        return True
    if not ishomog(signature[1:]) or bidegree(signature[0]) != bidegree(signature[1]):
        return False
    return True 

def deriv(poly):
    '''
    given a polynomial, returns a polynomial representing its image under the differential on the lambda algebra.
    '''
    if len(poly) == 0:
        return []
    if len(poly) == 1:
        if len(poly[0]) == 1:
            prod = []
            for i in range(1, poly[0][0] + 1):
                if binom(poly[0][0] - i, i) % 2 == 0:
                    continue
                prod.append((poly[0][0] - i, i - 1))
            return admiss(prod)
        return admiss(mul(deriv([(poly[0][0],)]), [poly[0][1:]]) + mul([(poly[0][0],)], deriv([poly[0][1:]]))) #TODO: change tuples to lists

    prod = []
    for mon in poly:
        prod += deriv([mon])
    
    return reduce(prod) 

#############################
# Helper functions: parsing #
#############################

def tokenize(string): #TODO: remove 'tuple(...)'; just make them lists
    '''
    Converts a string representation of a polynomial (e.g '1 0 + 2 1') to a list representation 
    (e.g [[1, 0], [2, 1]])
    '''
    prod = []
    thing = []
    num = ''

    if string.lower() == 'zero' or len(string) == 0:
        return prod

    for i in range(len(string)):
        if string[i] != ' ' and string[i] != '+':
            num += string[i]
            continue
        if string[i] == ' ':
            if string[i - 1] == '+':
                continue
            thing.append(int(num))
            num = ''
        if string[i + 1] == '+':
            prod.append(tuple(thing))
            thing = []
    
    thing.append(int(num))
    prod.append(tuple(thing))

    return prod

def detokenize(lst):
    '''
    Converts a list representation of a polynomial (e.g [[1, 0], [2, 1]]) to a string representation 
    (e.g '1 0 + 2 1')
    '''
    if len(lst) == 0:
        return 'Zero'

    mon = ''
    for i in range(len(lst)):
        for j in range(len(lst[i])):
            mon += f"{lst[i][j]}"
            if j != len(lst[i]) - 1:
                mon += ' '
        if i != len(lst) - 1:
            mon += ' + '
    return mon

########################
# The Polynomial class #
########################
class Polynomial:
    '''
    Class for polynomials in the lambda algebra.
    '''
    def __init__(self, signature):
        '''
        initializes polynomial object. signature is a list/tuple of list/tuples, 
        where each element is a tuple of numbers representing a monomial in the lambda algebra.
        '''
        self.sign = signature
        self.ishomogeneous = False
        self.bidegree = None

    def __str__(self):
        # if len(self.sign) == 0:
        #     return '0'
        # mon = ''
        # for i in range(len(self.sign)):
        #     for j in range(len(self.sign[i])):
        #         mon += f"{self.sign[i][j]}"
        #         if j != len(self.sign[i]) - 1:
        #             mon += ' '
        #     if i != len(self.sign) - 1:
        #         mon += ' + '
        # return mon
        return detokenize(self.sign)

    def weights(self):
        prod = []
        for i in self.sign:
            num = 0
            for j in i:
                num += math.ceil(j / 2)
            prod.append(num)
        return prod

    def __mul__(self, other):
        return Polynomial(mul(self.sign, other.sign))

    def __add__(self, other):
        return Polynomial(self.sign + other.sign)

    def admiss(self):
        return Polynomial(admiss(self.sign))

    def deriv(self):
        return Polynomial(deriv(self.sign))

    def __len__(self):
        return len(self.sign)

    def __lt__(self, other):
        return lexicographicorder2(self.sign, other.sign)        

    def __gt__(self, other):
        return not self < other and not self == other

    def __eq__(self, other):
        return self.sign == other.sign 

    def __le__(self, other):
        return not self > other

    def __ge__(self, other):
        return not self < other

######################################
# Helper functions: table generation #
######################################

def completecycle(poly, CTable):
    '''
    Given a polynomial in the lambda algebra, and dictionaries storing tagging information, 
    updates tagging information by trying to completing polynomial to a cycle; if the polynomial 
    is completed to a cycle, no updates are made, and otherwise updates input dictionaries accordingly.

    poly: polynomial in the lambda algebra (represented in list form)
    tags: dictionary containing information which tells us what monomials are tagging; tags[x] = y iff x tags y
    tagged: dictionary containing information which tells us what monomials are tagged: tagged[y] = x iff x tags y
    '''

    boundary = deriv(poly)

    if len(boundary) == 0: #if the polynomial is a cycle, then done
        return

    leadingterm = boundary[0]
    ltstr = detokenize([leadingterm])

    # if ltstr in CTable.tags:
    #     return 

    fragment = str(leadingterm[-1])
    
    if ltstr in CTable.table[leadingterm[0]][sum(leadingterm) - leadingterm[0]] and ltstr not in CTable.tagged:
        tagger = detokenize([poly[0]])
        taggee = ltstr
        CTable.tagged[taggee] = tagger
        CTable.tags[tagger] = taggee 
        if len(poly[0]) == 1 and leadingterm[-1] == 0:
            tagger = detokenize([poly[0] + (0,)])
            taggee = ltstr + ' 0'
            CTable.tagged[taggee] = tagger
            CTable.tags[tagger] = taggee 
        return

    for i in range(0,len(leadingterm)):
        #check to see if portion of leading term is in table
        portion = leadingterm[i:]
        portionstr = detokenize([portion])
        stem = sum(portion)
        if portionstr in CTable.table[portion[0]][stem - portion[0]]: #If CTable[i][j] were sets, this would go faster.
            if portionstr in CTable.tagged:
                preimage = [leadingterm[:i] + tokenize(CTable.tagged[portionstr])[0]]
                return completecycle(poly + preimage, CTable)

def longestmon(cell):
    '''
    Given a cell in the Curtis table, determines the length of the longest monomial in that cell.

    cell: list, (to be set) containing strings representing monomials in the lambda algebra.
    '''
    maxi = 0
    for i in cell:
        if len(i) > maxi:
            maxi = len(i)
    
    return maxi

def linscale(factor = 1):
    '''
    Returns the function "multiplication by (factor)".
    '''
    def prod(number):
        return factor*number
    return prod

def logscale(factor = 1):
    '''
    Returns the function "x \mapsto log (factor * x)"
    '''
    def prod(number):
        return math.log(factor * number)
    return prod

def makegrid(CTable, wfunc = linscale(1), hfunc = linscale(1)): #TODO: change wunit, hunit to functions wfunc, hfunc: \bbN \to \bbN instead, for a more sophisticated rule for scaling.
    '''
    Makes an empty grid in matplotlib that can house Curtis table.

    CTable: CurtisTable object.
    wfunc: function R -> R; takes in a number and outputs horizontal scaling factor associated to that number. Recommend it to be increasing.
    hfunc: function R -> R; takes in a number and outputs vertical scaling factor associated to that number. Recommend it to be increasing.
    '''
    width = sum([wfunc(i) for i in CTable.maxcellwidths])
    height = sum([hfunc(i) for i in CTable.maxcellsizes])

    #draw horizontal lines
    fig, ax = plt.subplots()
    # plt.rcParams['text.usetex'] = True

    ax.set_xlim([0, width])
    ax.set_ylim([-height,0])

    currheight = 0
    ax.axline((0, 0), (width, 0), linewidth = 0.5)
    ytix = []

    for i in range(len(CTable.maxcellsizes)):
        currheight += hfunc(CTable.maxcellsizes[i])
        ax.axline((0, -currheight), (width, -currheight), linewidth = 0.5)
        ytix.append(-(currheight - hfunc(CTable.maxcellsizes[i]) / 2))


    #draw vertical lines
    currwidth = 0
    ax.axline((0, 0), (0, height), linewidth = 0.5)
    xtix = []

    for i in range(len(CTable.maxcellwidths)):
        currwidth += wfunc(CTable.maxcellwidths[i])
        ax.axline((currwidth, 0), (currwidth, height), linewidth = 0.5)
        xtix.append(currwidth - wfunc(CTable.maxcellwidths[i]) / 2)

    xlabels = []
    ylabels = []

    for i in range(len(xtix)):
        xlabels.append(i)
        ylabels.append((2*i + 1, i + 1))

    ax.set_xticks(xtix, xlabels)
    ax.set_yticks(ytix, ylabels, rotation = 45)

    ax.xaxis.tick_top()
    ax.xaxis.set_label_position('top')
    ax.set_xlabel(r'$t - s$', labelpad = 8, color="blue")
    ax.set_ylabel(r'($2n - 1$, $n$)', labelpad = 10, rotation = 45, color="blue")
    
    return fig, ax

def gatherpositions(CTable, wfunc, hfunc):
    '''
    Gathers the position coordinates of each dot that will be placed in the Curtis table.

    CTable: CurtisTable object.
    wfunc: function R -> R; takes in a number and outputs horizontal scaling factor associated to that number. Recommend it to be increasing.
    hfunc: function R -> R; takes in a number and outputs vertical scaling factor associated to that number. Recommend it to be increasing.
    '''
    #gather the y-coordinates of the top of each row. The last entry will be the bottom of the bottommost row.
    lasttop = 0
    tops = [lasttop]
    for i in range(len(CTable.maxcellsizes)):
        lasttop -= hfunc(CTable.maxcellsizes[i])
        tops.append(lasttop)

    #gather the x-coordinates of the left of each column. The last entry will be the right of the rightmost row.
    lastleft = 0
    lefts = [lastleft]
    for i in range(len(CTable.maxcellwidths)):
        lastleft += wfunc(CTable.maxcellwidths[i])
        lefts.append(lastleft)

    xcoors = []
    ycoors = []
    positions = {}
    monomialat = {}

    for i in range(len(CTable.table)):
        for j in range(len(CTable.table[i])):
            top = tops[i]
            bottom = tops[i + 1]
            left = lefts[i + j]
            right = lefts[i + j + 1]

            xcoor = left + (right - left)/6
            yinc = (bottom - top)/(len(CTable.table[i][j]) + 1)
            ycoor = top + yinc

            for elt in CTable.table[i][j]:
                xcoors.append(xcoor)
                ycoors.append(ycoor)
                positions[elt] = (xcoor, ycoor)
                monomialat[(xcoor, ycoor)] = elt

                ycoor += yinc

    return xcoors, ycoors, positions, monomialat

def loadCTable(stem):
    '''
    If exists in directory, loads classical Curtis table, completed out to specified stem.

    stem: Int, specifying stem.
    '''
    if os.path.exists(f"ClassicalTables/ClassicalTableStem{stem}.py"):
        module = importlib.import_module(f"ClassicalTables.ClassicalTableStem{stem}")
        # module = __import__("ClassicalTables", fromlist=[f"ClassicalTableStem{stem}"])
        # module = module.ClassicalTableStem{stem}
        return CurtisTable(module.table, module.tags, module.tagged)

def gathercoordinatesE2(CTable):
    '''
    Gather the coordinates of afjdhjfkahs
    '''
    bidegrees = {}
    withbidegree = {}

    #gather the monomials which do not kill/are not killed, excluding the ones in the column furthest out.
    for i in range(len(CTable.table) - 1):
        for j in range(len(CTable.table[i]) - 1):
            for elt in CTable.table[i][j]:
                if elt not in CTable.tagged and elt not in CTable.tags and elt != '':
                    bd = bidegree(tokenize(elt)[0])
                    bidegrees[elt] = bd
                    if bd in withbidegree:
                        withbidegree[bd].append(elt)
                    else:
                        withbidegree[bd] = [elt]
    
    xcoords = []
    ycoords = []
    positions = {}
    monomialat = {}

    step = 0.2

    # for bd in withbidegree:
    #     lst = withbidegree[bd]
    #     count = len(lst)
    #     if count % 2 == 1:
    #         start = bd[0] - (count // 2)*step
    #     else:
    #         start = bd[0] - (count / 2 - 1)*step - step/2
    #     for i in range(count):
    #         xcoords.append(start)
    #         ycoords.append()

    for elt in bidegrees:
        if elt != '':
            ycoords.append(bidegrees[elt][1])
            if len(withbidegree[bidegrees[elt]]) % 2 == 1:
                start = bidegrees[elt][0] - (len(withbidegree[bidegrees[elt]]) // 2)*step
            else:
                start = bidegrees[elt][0] - (len(withbidegree[bidegrees[elt]]) / 2 - 1)*(step) - step/2
            for i in range(len(withbidegree[bidegrees[elt]])):
                if (start, bidegrees[elt][1]) not in monomialat:
                    xcoords.append(start)
                    positions[elt] = (start, bidegrees[elt][1])
                    monomialat[(start, bidegrees[elt][1])] = elt
                    break
                else:
                    start += step

    return xcoords, ycoords, positions, monomialat


#########################
# The CurtisTable class #
#########################

class CurtisTable: #TODO: Modify Curtis tables so that each cell is a set, not a list. (faster checking membership.)
    '''
    Class for Curtis tables.
    '''
    def __init__(self, table=[[['0', '']]], tags={}, tagged = {}):
        '''
        Initializer. Creates CurtisTable object, with attributes 'table' and 'arrows'.
        
        table: list of list of strings, storing entries in Curtis table. table[i] corresponds to row i, and table[i][j]
        is the i + jth box of row i. An element of table[i][j] is a monomial in the lambda algebra (represented as a list).

        tags: dictionary, with keys and values both being string representations of monomials. tags[x] = y means
        'a polynomial with leading term x has boundary with leading term y;' i.e 'x tags y'. In other words, if 
        x is a key in tags, then there is an arrow from x to y in the Curtis table.

        tags: dictionary, with keys and values both being string representations of monomials. tagged[x] = y means
        'a polynomial with leading term y has boundary with leading term x;' i.e 'y tags x'. In other words, if 
        y is a key in tagged, then there is an arrow from x to y in the Curtis table.
        '''
        for i in table:
            if not isinstance(i, list):
                raise LambdaEvaluationError('table malformed.')
        self.table = table
        self.stem = len(table) - 1
        self.maxcellsizes = [max([len(x) for x in table[i]]) for i in range(len(table))]
        self.maxcellwidths = [max([longestmon(table[i - j][j]) for j in range(i + 1)]) for i in range(self.stem + 1)]
        self.tags = tags
        self.tagged = tagged

    def save(self):
        '''
        Saves data in Curtis table to Python file, in addition to a png of the current table. If data for table of the same stem already 
        exists in directory, does nothing. Data can then be loaded into CurtisTable object using function 'load(stem)'.

        Filename of file where data is stored is "ClassicalTableStem(the stem the table is completed until).py"
        '''
        if os.path.exists(f"ClassicalTables/ClassicalTableStem{self.stem}.py"):
            return
        prod = open(f"ClassicalTables/ClassicalTableStem{self.stem}.py", "w")
        prod.write(f"table = {x.table}\n")
        prod.write(f"tags = {x.tags}\n")
        prod.write(f"tagged = {x.tagged}\n")
        prod.close
        return

    def saveimage(self, show_survivors = True, show_killed = False, dpi=500):
        if os.path.exists(f"ClassicalTablesImages/ClassicalTableStem{self.stem}, ss={show_survivors}, sk={show_killed}, dpi={dpi}.png"):
            return
        fig, _ = self.display(output = False, show_survivors=show_survivors, show_killed = show_killed)
        fig.savefig(f"ClassicalTablesImages/ClassicalTableStem{self.stem} ss={show_survivors}, sk={show_killed}, dpi={dpi}.png", 
            format = 'png', 
            dpi = dpi, 
            bbox_inches = 'tight'
            )
    
    def expand(self):
        '''
        Expands classical Curtis table by one.
        '''
        self.stem += 1
        self.table[0].append([])
        self.table.append([]) #add new row
        
        if self.stem == 1: #may need to implement special case for stem = 1
            pass
        
        #fill out column self.stem, aside from top 2 cells.
        for i in range(self.stem, 1, -1): #add cell onto row i (starting at bottom)
            prod = []
            #want to look at all permanent cycles living in stem (self.stem - i)...
            #of which there are i + 1
            for j in range(0, min(2*i + 1, self.stem - i + 1)): #?? high chance of bugs
                for mon in self.table[j][self.stem - i - j]: #バグ可能性大 -> remember that you only want to check stuff in \Lambda(2n - 1)!
                    if (mon not in self.tagged and mon not in self.tags) or (mon in self.tagged and tokenize(self.tagged[mon])[0][0] >= 2*i + 1):
                        if mon != '':
                            prod.append(str(i) + ' ' + mon)
                        else:
                            prod.append(str(i))
                # except Exception as e:
                #     print(e)
                #     print((i, j))
            self.table[i].append(prod)

        #compute differentials in new column, top down.
        for i in range(2, self.stem + 1):
            for mon in self.table[i][-1]:
                # if mon == '2 3 1':
                #     print('here')
                completecycle(tokenize(mon), self)

        #TODO: fill out cell in row 2
        prod = []
        if self.stem > 1:
            for mon in self.table[1][self.stem - 2]:
                if (mon not in self.tagged and mon not in self.tags) or (mon in self.tagged and tokenize(self.tagged[mon])[0][0] >=  3):
                    prod.append(str(1) + ' ' + mon)
            for mon in self.table[2][self.stem - 3]:
                if mon not in self.tags:
                    prod.append(str(1) + ' ' + mon)
        else:
            prod.append(str(1))
            prod.append(str(1) + ' 0')
        
        self.table[1].append(prod)
        
        #update maxcellsizes
        for i in range(len(self.maxcellsizes)):
            if len(self.table[i][-1]) > self.maxcellsizes[i]:
                self.maxcellsizes[i] = len(self.table[i][-1])

        self.maxcellsizes.append(2) #biggest size in the newly made row is always 2
        self.maxcellwidths.append(max([longestmon(self.table[j][self.stem - j]) for j in range(self.stem + 1)])) 

    def __str__(self):
        '''
        str method. Helps to print the table and tags in a semi-readable fashion. Useful for debugging.
        '''
        prod = f'stem {self.stem}\ntable:\n'
        for i in range(len(self.table)):
            prod += f'{self.table[i]}\n'
        prod += '\ntags:\n'
        for i in self.tags:
            prod += f'{(i, self.tags[i])}, '
        return prod

    def display(self, output = True, wfunc = logscale(2), hfunc = linscale(), show_survivors = True, show_killed = False, hover_to_show = False):
        '''
        Returns and displays Curtis table as interactive matplotlib figure. Returns: (fig, ax), matplotlib figure and axis objects
        
        output: boolean value; specifies whether or not to display the table. default is True.
        wfunc: function R -> R; takes in a number and outputs horizontal scaling factor associated to that number. Recommend it to be increasing.
        hfunc: function R -> R; takes in a number and outputs vertical scaling factor associated to that number. Recommend it to be increasing.
        show_survivors: boolean value; specifies whether or not to display the names of elements which are not tags/are not tagged. (permanent cycles)
        show_killed: boolean value; specifies whether of not to display the names of elements which are tags/get tagged.
        hover_to_show: boolean value; specifies whether or not to enable hovering to show names of elements. If enabled, elements will only show their names
        when hovered over with mouse.
        '''
        fig, ax = makegrid(self, wfunc, hfunc)
        xcoors, ycoors, positions, monomialat = gatherpositions(self, wfunc, hfunc)
        #change the colors of those that aren't tagged
        colors = []
        for i in range(len(xcoors)):
            if monomialat[(xcoors[i], ycoors[i])] not in self.tags and monomialat[(xcoors[i], ycoors[i])] not in self.tagged:
                colors.append('indianred')
            else:
                colors.append('darkgray')

        scatter = plt.scatter(xcoors, ycoors, s = 1.5, color = colors)

        #add tags to each dot, telling what monomial it is.
        width = sum([wfunc(self.maxcellwidths[i]) for i in range(len(self.maxcellwidths))])
        height = sum([hfunc(self.maxcellsizes[i]) for i in range(len(self.maxcellsizes))])

        for i in self.tags:
            try:
                start = positions[i]
                end = positions[self.tags[i]]

                # dx, dy = end[0] - start[0], end[1] - start[1]
                ax.annotate(
                    "", xy=end, xytext=start,
                    arrowprops=dict(arrowstyle="->", lw=0.7, color='darkgray')
                    )
            except:
                continue
        
        if hover_to_show == False:
            for i in range(len(xcoors)):
                if (monomialat[(xcoors[i],ycoors[i])] in self.tags or monomialat[(xcoors[i],ycoors[i])] in self.tagged) and show_killed == False:
                    continue
                if monomialat[(xcoors[i],ycoors[i])] not in self.tags and monomialat[(xcoors[i],ycoors[i])] not in self.tagged and show_survivors == False:
                    continue
                annot = ax.annotate(
                    text = monomialat[(xcoors[i], ycoors[i])], 
                    xy = (xcoors[i], ycoors[i]),
                    xytext = (width/20, -height/20), #TODO: find less arbitrary way of setting this.
                    textcoords='offset points',
                    color = colors[i],
                    fontsize = 'small'
                    )
        elif hover_to_show == True:
            annot = ax.annotate(
                text = '',
                xy = (0,0),
                xytext = (width/20, -height/20), #TODO: find less arbitrary way of setting this.
                textcoords='offset points',
                weight = 'bold',
                fontsize = 'small'
                )
                
            annot.set_visible(False)

            def hover(event):
                is_visible = annot.get_visible()
                if event.inaxes == ax:
                    is_containing, index = scatter.contains(event)
                    if is_containing:
                        loc = scatter.get_offsets()[index['ind'][0]]
                        if show_killed == False and (monomialat[tuple(loc)] in self.tags or monomialat[tuple(loc)] in self.tagged):
                            pass
                        elif show_survivors == False and monomialat[tuple(loc)] not in self.tags and monomialat[tuple(loc)] not in self.tagged:
                            pass
                        annot.xy = loc
                        annot.set_text(monomialat[tuple(loc)])
                        annot.set(color = colors[index['ind'][0]])
                        annot.set_visible(True)

                        fig.canvas.draw_idle()
                    else:
                        if is_visible:
                            annot.set_visible(False)
                            fig.canvas.draw_idle()

            fig.canvas.mpl_connect('motion_notify_event', hover)
        
        if output:
            plt.show()

        return fig, ax

    def displayE2(self, output = True, show_names = True, hover_to_show = False):
        '''
        Uses the data in the Curtis table to display the E2 page of the classical Adams spectral sequence.
        '''
        #gather data needed to determine grid size.

        xcoords, ycoords, positions, monomialat = gathercoordinatesE2(self)

        width = self.stem
        height = max(ycoords) + 1

        #generate grid
        fig, ax = plt.subplots()

        ax.set_xlim([-0.5, width])
        ax.set_ylim([-0.5, height])

        ax.set_xlabel(r'$t - s$', labelpad = 8, color="blue")
        ax.set_ylabel(r'$s$', labelpad = 10, rotation = 0, color="blue")

        #draw vertical lines

        start = 0
        ax.axline((start,0), (start, height), linewidth = 1, color = 'black', zorder = -1)
        for i in range(width):
            start += 1
            ax.axline((start,0), (start, height), linewidth = 0.5, color = 'lightgray', zorder = -1)
            

        #draw horizontal lines
        start = 0
        ax.axline((0, start), (width, start), linewidth = 1, color = 'black', zorder = -1)
        for i in range(height + 1):
            start += 1
            ax.axline((0, start), (width, start), linewidth = 0.5, color = 'lightgray', zorder = -1)
        
        scatter = plt.scatter(xcoords, ycoords, s= 5, color = 'black', zorder = 1)

        if show_names:
            if hover_to_show == False:
                for i in range(len(xcoords)):
                    annot = ax.annotate(
                        text = monomialat[(xcoords[i], ycoords[i])], 
                        xy = (xcoords[i], ycoords[i]),
                        xytext = (0, height/20), #TODO: find less arbitrary way of setting this.
                        textcoords='offset points',
                        color = 'black',
                        fontsize = 'xx-small',
                        ha = 'center'
                        )
            elif hover_to_show == True:
                annot = ax.annotate(
                    text = '',
                    xy = (0,0),
                    xytext = (0, height/5), #TODO: find less arbitrary way of setting this.
                    textcoords='offset points',
                    weight = 'bold',
                    fontsize = 'small',
                    ha = 'center'
                    )
                    
                annot.set_visible(False)

                def hover(event):
                    is_visible = annot.get_visible()
                    if event.inaxes == ax:
                        is_containing, index = scatter.contains(event)
                        if is_containing:
                            loc = scatter.get_offsets()[index['ind'][0]]
                            annot.xy = loc
                            annot.set_text(monomialat[tuple(loc)])
                            annot.set(color = 'black')
                            annot.set_visible(True)

                            fig.canvas.draw_idle()
                        else:
                            if is_visible:
                                annot.set_visible(False)
                                fig.canvas.draw_idle()

                fig.canvas.mpl_connect('motion_notify_event', hover)

        if output:
            plt.show()

        return fig, ax

if __name__ == '__main__':
    doctest.testmod()
    import test
    import time

    # x = 'fitness'
    # print({x})
    # x = loadCTable(26)
    # x.display(show_killed=True, hover_to_show=True)
    # n = 32

    # x = CurtisTable()

    # start = time.perf_counter()
    # for i in range(n - x.stem):
    #     x.expand()
    #     end = time.perf_counter()
    #     print(f"stem {x.stem - 1} to stem {x.stem}: {end - start:0.4f} seconds")
    #     x.save()
    #     start = time.perf_counter()
    
    # x.save()
    x = loadCTable(26)

    # x.display(show_survivors = False)
    x.displayE2(hover_to_show=True)

    replon = False

    while replon:
        print("What do you want to do?")
        x = input('in> ')
        if x.lower() == 'quit':
            break
        if x.lower() == 'simplify':
            flag = True
            while flag:
                print("Input polynomial you want to simplify.")
                y = input('in> ')
                try:
                    poly = tokenize(y)
                    given = Polynomial(poly)
                    print(f"\n {given} \n \n in admissible form is \n \n {given.admiss()} \n")
                except:
                    try: 
                        if y.lower() == 'quit':
                            flag = False
                    except:
                        raise LambdaSyntaxError('Input polynomial not in correct form.')
        elif x.lower() == 'compute differential':
            flag = True
            while flag:
                print("Input polynomial you want to compute the differential of.")
                y = input('in> ')
                try:
                    poly = tokenize(y)
                    given = Polynomial(poly)
                    boundary = given.deriv()
                    print(f"\n The differential of \n \n {given}   (motivic weights: {given.weights()}) \n \n is \n \n {boundary}   (motivic weights: {boundary.weights()}). \n ")
                except:
                    try:
                        if y.lower() == 'quit':
                            flag = False
                    except:
                        raise LambdaSyntaxError('Input polynomial not in correct form.')
        else:
            raise LambdaEvaluationError('Command not defined.')