#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 26 17:29:18 2022

@author: Roland Hansson
"""

import numpy as np
import sympy as sym
import matplotlib.pyplot as plt
import matplotlib.colors as col
import sys
from NewtonFractal2D import fractal2D


import unittest

class TestIdentity(unittest.TestCase):

#############################
### TESTS NEWTON'S METHOD ###
#############################
    

################ NOTE: NEEDS TO CHECK ANSWERS ARE CORRECT #############


#### FIRST SET OF EQUATIONS
    def test_newton1(self):
        a, b = sym.Symbol('a'), sym.Symbol('b')        
        eq1, eq2 = a**3 - 3*a*b**2 - 1, 3*a**2*b - b**3  
        p = fractal2D(eq1, eq2, 100)  
        
        
        result = p.newtonsMethod( [1,0] )          #Directly inputting a root
        expected = [1, 0]
        rounded = [np.round(x,5) for x in result]
        self.assertEqual(rounded[0] , expected[0]) #assertAlmostEqual doesn't work with arrays
        self.assertEqual(rounded[1] , expected[1]) #because apparently arrays don't have rounding
        
        
        result = p.newtonsMethod( [-2,-3] )        
        expected = [-0.5, -0.86603]             #Testing negative input values
        rounded = [np.round(x,5) for x in result]
        self.assertEqual(rounded[0] , expected[0]) 
        self.assertEqual(rounded[1] , expected[1])        
        

        p = fractal2D(eq1, eq2, 100)        
        result = p.newtonsMethod( [20,30] )        
        expected = [-0.5, 0.86603]
        rounded = [np.round(x,5) for x in result]
        self.assertEqual(rounded[0] , expected[0]) 
        self.assertEqual(rounded[1] , expected[1])        
            
        
        
        
#### SECOND SET OF EQUATIONS        
    def test_newton2(self):
        a, b = sym.Symbol('a'), sym.Symbol('b')       
        eq1, eq2 = a**3 - 3*a*b**2 - 2*a - 2, 3*a**2*b - b**3 - 2*b  
        p = fractal2D(eq1, eq2, 100)  
        
        
        result = p.newtonsMethod( [1.76929,0] )    #Directly inputting a root    
        expected = [1.76929, 0]
        rounded = [np.round(x,5) for x in result]
        self.assertEqual(rounded[0] , expected[0]) #assertAlmostEqual doesn't work with arrays
        self.assertEqual(rounded[1] , expected[1]) #because apparently arrays don't have rounding
        
        
        result = p.newtonsMethod( [-2,-3] )        #Testing negative input values
        expected = [-0.88465, -0.58974]
        rounded = [np.round(x,5) for x in result]
        self.assertEqual(rounded[0] , expected[0]) 
        self.assertEqual(rounded[1] , expected[1])        
        

        p = fractal2D(eq1, eq2, 100)        
        result = p.newtonsMethod( [20,30] )        
        expected = [-0.88465, 0.58974]
        rounded = [np.round(x,5) for x in result]
        self.assertEqual(rounded[0] , expected[0]) 
        self.assertEqual(rounded[1] , expected[1])                
      
        
      
        
#### THIRD SET OF EQUATIONS      
    def test_newton3(self):
        a, b = sym.Symbol('a'), sym.Symbol('b')        
        eq1 = a**8 - 28*a**6*b**2 + 70*a**4*b**4 + 15*a**4 -28*a**2*b**6 -90*a**2*b**2 +b**8 -16
        eq2 = 8*a**7*b - 56*a**5*b**3 + 56*a**3*b**5 - 60*a**3*b - 8*a*b**7 - 60*a*b**3  
        p = fractal2D(eq1, eq2, 100)  
        
        
        result = p.newtonsMethod( [-1,0] )         #Directly inputting a root
        expected = [-1, 0]
        rounded = [np.round(x,5) for x in result]
        self.assertEqual(rounded[0] , expected[0]) #assertAlmostEqual doesn't work with arrays
        self.assertEqual(rounded[1] , expected[1]) #because apparently arrays don't have rounding
        
        
        result = p.newtonsMethod( [-2,-3] )        #Testing negative input values
        expected = [-1.49685, -1.91762]
        rounded = [np.round(x,5) for x in result]
        self.assertEqual(rounded[0] , expected[0]) 
        self.assertEqual(rounded[1] , expected[1])        
        

        p = fractal2D(eq1, eq2, 100)        
        result = p.newtonsMethod( [20,30] )        
        expected = [1.49685, 1.91762]
        rounded = [np.round(x,5) for x in result]
        self.assertEqual(rounded[0] , expected[0]) 
        self.assertEqual(rounded[1] , expected[1])      




#########################
### TESTS GRAPH MAKER ###
#########################        



        
    def test_shade_colour(self):
        a, b = sym.Symbol('a'), sym.Symbol('b') 
        eq1, eq2 = a**3 - 3*a*b**2 - 1, 3*a**2*b - b**3
        f = fractal2D(eq1, eq2, 100)        
        f.A = np.array([[0,1,2,3],[0,1,2,3],[0,1,2,3]])
        f.B = np.array([[5,5,5,5],[6,6,6,6],[6,7,8,9]])
        f.nrIterations = 9
        
        
        #### SHADE COLOUR
        base = (2.,2.,2.)
        result = f.shade_color( base, 4, 4, 9) 
        expected = (2.,2.,2.)                   #Unchanged if equal to minimum
        self.assertEqual( result , expected ) 
               
        base = (1.,1.,1.)
        result = f.shade_color( base, 10, 6, 10)
        expected = (0.2,0.2,0.2)
        self.assertEqual( result , expected ) 
        
        
        
        
    def test_make_cmap(self):
        a, b = sym.Symbol('a'), sym.Symbol('b') 
        eq1, eq2 = a**3 - 3*a*b**2 - 1, 3*a**2*b - b**3
        f = fractal2D(eq1, eq2, 100)        
        f.A = np.array([[0,1,2,3],[0,1,2,3],[0,1,2,3]])
        f.B = np.array([[5,5,5,5],[6,6,6,6],[6,7,8,9]])
        f.nrIterations = 9
        
        
        #### MAKE CMAP
        cmap = f.make_cmap()
        result = cmap(0.)
        #print(result)
        expected = (1.0, 0.5, 0.5, 1.0)
        self.assertEqual( result , expected ) 
        
        
        
        
    def test_make_graph(self):
        print("MAKE GRAPH")
        a, b = sym.Symbol('a'), sym.Symbol('b') 
        eq1, eq2 = a**3 - 3*a*b**2 - 1, 3*a**2*b - b**3
        f = fractal2D(eq1, eq2, 100)   
        
        #Makes a 3x4 grid, red all one color, green blue and yellow getting darker. 
        f.A = np.array([[0,1,2,3],[0,1,2,3],[0,1,2,3]])
        f.B = np.array([[5,5,5,5],[5,6,6,6],[5,7,8,9]])
        f.graph()
        
        #Makes a 3x4 grid, red all one color, green blue and yellow getting darker. 
        f.A = np.array([[0.,1.,2.,3.],[0.,1.,2.,3.],[0.,1.,2.,3.]])
        f.B = np.array([[5,5,5,5],[5,6,6,6],[5,7,8,9]])
        f.graph()
        
            
        
if __name__=='__main__':
    unittest.main()



# Testing: Init

# Testing: NewtonMethod

# Testing: Root Finding

# Testing: Graph


