module Tests.Cholesky

import Cholesky

import Data.Vect


assertEq : Eq a => (given : a) -> (expected : a) -> IO ()
assertEq g e = if g == e
    then putStrLn "Test Passed"
    else putStrLn "Test Failed"

assertNotEq : Eq a => (given : a) -> (expected : a) -> IO ()
assertNotEq g e = if not (g == e)
    then putStrLn "Test Passed"
    else putStrLn "Test Failed"


assertTrue : Bool -> IO ()
assertTrue True = putStrLn "Test Passed"
assertTrue False = putStrLn "Test Failed"

 
test3By3Cholesky : IO ()
test3By3Cholesky = assertEq expected calculated
	where
		expected : Matrix 3 3 Double
		expected = [[5.0, 0.0, 0.0], [3.0, 3.0, 0.0], [-1.0, 1.0, 3.0]]

		calculated : Matrix 3 3 Double
		calculated = cholesky [[25.0, 15.0, -5.0], [15.0, 18.0, 0.0], [-5.0, 0.0, 11.0]] 


test4By4Cholesky : IO ()
test4By4Cholesky = do 
		assertTrue $ isEqual expected calculated 0.0001
	where
		expected : Matrix 4 4 Double
		expected = [[4.242640687119285, 0, 0, 0], [5.185449728701349, 6.565905201197403, 0, 0], 
								[12.72792206135786, 3.046038495400855, 1.64974224790907, 0], 
								[9.899494936611665, 1.624553864213789, 1.849711005231382, 1.392621247645587]]

		calculated : Matrix 4 4 Double
		calculated = cholesky [[18.0, 22.0, 54.0, 42.0], [22.0, 70.0, 86.0, 62.0],
                						[54.0, 86.0, 174.0, 134.0], [42.0, 62.0, 134.0, 106.0]] 

