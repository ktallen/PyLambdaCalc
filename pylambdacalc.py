import doctest
import scipy.special

#####################
# Custom Exceptions #
#####################
class LambdaError(Exception):
    pass

class LambdaSyntaxError(LambdaError):
    pass

class LambdaEvaluationError(LambdaError):
    pass

####################
# Helper functions #
####################
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
        return admiss(mul(deriv([(poly[0][0],)]), [poly[0][1:]]) + mul([(poly[0][0],)], deriv([poly[0][1:]])))

    prod = []
    for mon in poly:
        prod += deriv([mon])
    
    return reduce(prod) #maybe should call admiss

def tokenize(string):
    prod = []
    thing = []
    num = ''

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

def str_to_poly(string):
    '''
    Converts string representation of polynomial to list/tuple representation of polynomial.
    '''
    if str_to_poly == '' or str_to_poly == '0':
        return []

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
        if len(self.sign) == 0:
            return '0'
        mon = ''
        for i in range(len(self.sign)):
            for j in range(len(self.sign[i])):
                mon += f"{self.sign[i][j]}"
                if j != len(self.sign[i]) - 1:
                    mon += ' '
            if i != len(self.sign) - 1:
                mon += ' + '
        return mon

    def __mul__(self, other):
        # sign = []
        # for mon1 in self.sign:
        #     for mon2 in other.sign:
        #         sign.append(mon1 + mon2)
        # return Polynomial(sign)
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

if __name__ == '__main__':
    doctest.testmod()
    import test

    replon = True

    while replon:
        print("What do you want to do?")
        x = input('in> ')
        if x.lower() == 'quit':
            break
        if x.lower() == 'simplify':
            print("Input polynomial you want to simplify.")
            y = input('in>')
            try:
                poly = tokenize(y)
                given = Polynomial(poly)
                print(f"\n {given} \n \n in admissible form is \n \n {given.admiss()} \n")
            except:
                raise LambdaSyntaxError('Input polynomial not in correct form.')
        elif x.lower() == 'compute differential':
            print("Input polynomial you want to compute the differential of.")
            y = input('in>')
            try:
                poly = tokenize(y)
                given = Polynomial(poly)
                print(f"\n The differential of \n \n {given} \n \n is \n \n {given.deriv()}. \n ")
            except:
                raise LambdaSyntaxError('Input polynomial not in correct form.')
        else:
            raise LambdaEvaluationError('Command not defined.')