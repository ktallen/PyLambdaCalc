import doctest
import scipy.special
import math

# from helpers import *

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

def admiss(poly): #think about the process of making a monomial admissible
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
    Returns bidegree of a monomial
    '''
    return len(monomial), sum(monomial) + len(monomial)

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
    
    return reduce(prod) #maybe should call admiss

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

    if string.lower() == 'zero':
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

######################################
# Helper functions: table generation #
######################################

def completecycle(poly, CTable): #not right; prob need to account for case when something is tagged by something in the same column
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

    if ltstr in CTable.tags:
        return 

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
            CTable.tags[tagger] = taggee #im kinda hardcoding this... and its wrong! see notes

    #check to see if 
    for i in range(1,len(leadingterm)):
        #check to see if portion of leading term is in table
        portion = leadingterm[i:]
        portionstr = detokenize([portion])
        stem = sum(portion)
        if portionstr in CTable.table[portion[0]][stem - portion[0]]:
            # if i == 0 and portionstr not in CTable.tagged:
            #     tagger = detokenize([poly[0]])
            #     taggee = portionstr
            #     CTable.tagged[taggee] = tagger
            #     CTable.tags[tagger] = taggee 
            #     return
            if portionstr in CTable.tagged:
                preimage = [leadingterm[:i] + tokenize(CTable.tagged[portionstr])[0]]
                return completecycle(poly + preimage, CTable)

    # if leadingterm[-1] == 0: #this is sus; might modify stuff i dont want it to.
    #     trunc = leadingterm[:-1]
    #     completecycle([trunc], CTable)
    #     truncstr = detokenize([trunc])
    #     if truncstr in CTable.tagged:
    #         CTable.tagged[ltstr] = CTable.tagged[truncstr] + ' 0'
    #         CTable.tags[CTable.tagged[ltstr]] = ltstr
        

    # for i in range(len(leadingterm) - 2, -1, -1): #maybe should be range(..., 0, -1)
    #     fragment = str(leadingterm[i]) + ' ' + fragment
    #     if fragment in CTable.tagged:
    #         try:
    #             preimage = tokenize(CTable.tagged[fragment])
    #             summand = [leadingterm[:i] + preimage[0]]
    #             return completecycle(poly + summand, CTable)
    #         except Exception as e:
    #             print(e)
    #             print(preimage)
    #             print(leadingterm[:i])
    
    # # fragment = str(leadingterm[0]) + ' ' + fragment
    # tagger = detokenize([poly[0]])
    # taggee = detokenize([leadingterm])

    # CTable.tagged[taggee] = tagger
    # CTable.tags[tagger] = taggee
    # return

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

class CurtisTable:
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
            if not isinstance(i, list) or len(i) != len(table):
                raise LambdaEvaluationError('table malformed.')
        self.table = table
        self.stem = 0
        self.maxcellsizes = [max([len(x) for x in table[i]]) for i in range(len(table))]
        self.tags = tags
        self.tagged = tagged
    
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
                try:
                    for mon in self.table[j][self.stem - i - j]: #バグ可能性大 -> remember that you only want to check stuff in \Lambda(2n - 1)!
                        if (mon not in self.tagged and mon not in self.tags) or (mon in self.tagged and tokenize(self.tagged[mon])[0][0] >= 2*i + 1):
                            if mon != '':
                                prod.append(str(i) + ' ' + mon)
                            else:
                                prod.append(str(i))
                except Exception as e:
                    print(e)
                    print((i, j))
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
                self.maxcellsizes[i] = self.maxcellsizes[i]

        self.maxcellsizes.append(2) #biggest size in the newly made row is always 2

    def __str__(self):
        prod = ''
        for i in range(len(self.table)):
            prod += f'{self.table[i]}\n'

        return prod
        

if __name__ == '__main__':
    doctest.testmod()
    import test

    # x = CurtisTable()
    # n = 13

    # for i in range(n):
    #     x.expand()
    
    # print(f"\nstem {n}:")
    # # print(f"{x.table}")
    # print(x)
    # print('tags:')
    # print([(j, x.tags[j]) for j in x.tags])
    # print('\n')

    # x = [1, 2, 3, 4]
    # print(x[:-1])

    replon = True

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