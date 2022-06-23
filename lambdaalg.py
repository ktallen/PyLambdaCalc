import doctest

# NO ADDITIONAL IMPORTS ALLOWED!
# You are welcome to modify the classes below, as well as to implement new
# classes and helper functions as necessary.

# def tokenize(string):
#     '''
#     Given a string, returns a list of strings. These strings should represent the list split into 
#     parentheses, numbers, and variables.
#     '''
#     prod = [] 
#     for i in string.split():
#         if '(' not in i and ')' not in i:
#             prod.append(i)
#         elif '(' in i:
#             if i.count('(') == 1:
#                 prod += [i[0], i[1:]]
#             else:
#                 prod += [i[0]]
#                 prod += tokenize(i[1:])
#         elif ')' in i:
#             if i.count(')') == 1:
#                 prod += [i[:-1], i[-1]]
#             else:
#                 prod += tokenize(i[:-1])
#                 prod += [i[-1]]
#     return prod

# def parse(token):
#     '''
#     Given a tokenized expession, parses it into a variable expression.
#     '''
#     def parse_expression(index):
#         if token[index] == '+' or token[index] == '-' or token[index] == '*' or token[index] == '/' or token[index] == ')':
#             index += 1
#         try:
#             thing = float(token[index])
#             return Num(thing), index + 1
#         except:
#             if token[index] != '(':
#                 return Var(token[index]), index + 1
#             else:
#                 left, next1 = parse_expression(index + 1)
#                 optype = token[next1]
#                 right, next2 = parse_expression(next1)
#                 if optype == '+':
#                     return Add(left, right), next2 + 1
#                 if optype == '-':
#                     return Sub(left, right), next2 + 1
#                 if optype == '*':
#                     return Mul(left, right), next2 + 1
#                 if optype == '/':
#                     return Div(left, right), next2 + 1

#     parsed_expression, _ = parse_expression(0)
#     return parsed_expression

# def sym(string):
#     return parse(tokenize(string))

########################
# Error messages
########################

def tokenize(string):
    pass

def binomcalc(n, m):
    '''
    Calculates binomial coefficients mod 2.
    '''
    if m == 0 and n == 0:
        return 1
    
    

class Polynomial:
    def __init__(self, signature):
        '''initializer. creates a monomial representing a given list of numbers; assumes
        is inadmissible.'''
        self.name = signature
        self.admissible = False
        self.expanded = False
        self.tag = None
        self.bidegree = None

    def __str__(self):
        '''
        returns string representation of monomial
        '''
        return self.name

    def __repr__(self):
        '''
        returns machine readable representation of variable.
        '''
        return 'Monomial(' + repr(self.name) + ')'
    
    def __lt__(self, other):
        pass

    def __gt__(self, other):
        pass

    def __eq__(self, other):
        pass

    def __ne__(self, other):
        pass

    def __le__(self, other):
        pass

    def __ge__(self, other):
        pass

class Monomial(Polynomial):
    '''
    Class for monomials in the lambda algebra.
    '''
    def __init__(self, signature):
        '''initializer. creates a monomial representing a given list of numbers; assumes
        is inadmissible.'''
        self.name = signature
        self.admissible = False

    def __str__(self):
        '''
        returns string representation of monomial
        '''
        return self.name

    def __repr__(self):
        '''
        returns machine readable representation of variable.
        '''
        return 'Monomial(' + repr(self.name) + ')'
    def isadmiss(lst):
        '''
        helper function; checks to see if a list of numbers representing a monomial in the
        lambda algebra is admissible.
        '''
        if not isinstance(lst, (list, tuple)):
            raise TypeError('Tried to check if nonlist is admissible')        

    def admiss(self):
        '''
        makes 
        '''
        if self.admissible == True:
            return
        for i in range(len(self.name)):
            pass
    
    def __lt__(self, other):
        pass

    def __gt__(self, other):
        pass

    def __eq__(self, other):
        pass

    def __ne__(self, other):
        pass

    def __le__(self, other):
        pass

    def __ge__(self, other):
        pass

class Generator(Monomial):
    def __init__(self, n):
        if not isinstance(n, int):
            raise TypeError("Tried to create generator with nonintegral degree.")
        self.index = n

    def deriv(self):
        if self.index <= 0:
            return Generator(-1) #Generator(-1) will be what I treat as 0.  

    def __lt__(self, other):
        if self.index < other.index:
            return True
        return False

    def __gt__(self, other):
        pass

    def __eq__(self, other):
        pass

    def __ne__(self, other):
        pass

    def __le__(self, other):
        pass

    def __ge__(self, other):
        pass
    

class BinOp(Polynomial):
    pass

class Add(BinOp):
    pass

class Mul(BinOp):
    pass



if __name__ == '__main__':
    doctest.testmod()
    import test
    # t = Mul(Div(Num(2), Num(2)), Var('x'))
    # print(t.simplify())

    # print(result)
    # print(' ')
    # print(expected)
    # print(sym('((z * 3) + 0)').deriv('z').simplify())