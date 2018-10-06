#!/usr/bin/env python 


class TestProperty(object):

    __hello__ = None
    @property
    def plus(self):

        return self.__hello__

    @plus.setter
    def plus(self, *args):

       #if value is True:
       #    self.__hello__ = True
       #else:
       #    self.__hello__ = args 
       self.__hello__ = args 


class Structure(object):
    s = 1.0


class loop(object):
    
    def lp(self, n, func):
       if n <5:
         l = True
       else: 
         l = False
       while l:
           n += 1
	 
           func.s = n
           print n, func.s
	   l=self.lp(n, func)
       return False 

a=loop()
b=Structure()
a.lp(1, b)
print b.s
#a= TestProperty()
#a.plus = True, 1,2,3,4
#a.plus = True  
#b=a.plus
#print ''.join(a.plus[0])
