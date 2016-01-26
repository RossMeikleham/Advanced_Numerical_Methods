module Cholesky
 
import Data.Vect
 
Matrix : Nat -> Nat -> Type -> Type
Matrix m n t = Vect m (Vect n t)
 
 
zeros : (m : Nat) -> (n : Nat) -> Matrix m n Double
zeros m n = replicate m (replicate n 0.0)
 
 
indexM : (Fin m, Fin n) -> Matrix m n t -> t
indexM (i, j) a = index j (index i a)
 
 
replaceAtM : (Fin m, Fin n) -> t -> Matrix m n t -> Matrix m n t
replaceAtM (i, j) e a = replaceAt i (replaceAt j e (index i a)) a
 

isEqualV : Vect m Double -> Vect m Double -> Double -> Bool
isEqualV v1 v2 delta = foldl (\eqs,(x, y) => eqs && abs(x - y) < delta) True (zip v1 v2)

isEqual : Matrix m n Double -> Matrix m n Double -> Double -> Bool
isEqual m1 m2 delta = foldl (\eqs,(v1, v2) => eqs && isEqualV v1 v2 delta) True (zip m1 m2) 

 
get : Matrix m m Double -> Matrix m m Double -> (Fin m, Fin m) -> Double
get a l (i, j) {m} =  if i == j then sqrt $ indexM (j, j) a - dot
             else if i > j  then (indexM (i, j) a - dot) / indexM (j, j) l 
             else 0.0
 
  where        
        -- Obtain indicies 0 to j -1 
        ks : List (Fin m)
        ks = case (findIndices (\_ => True) a) of
          [] => []
          (x::xs) => init (x::xs) 
 
        dot : Double
        dot = sum [(indexM (i, k) l) * (indexM (j, k) l) | k <- ks]
 
 
updateL : Matrix m m Double -> Matrix m m Double -> (Fin m, Fin m) -> Matrix m m Double
updateL a l idx = replaceAtM idx (get a l idx) l
 
 
cholesky : Matrix m m Double -> Matrix m m Double
cholesky a {m} = 
    foldl (\l',i => 
        foldl (\l'',j => updateL a l'' (i, j)) l' (js i)) 
          l is
  where  l = zeros m m
 
         is : List (Fin m) 
         is = findIndices (\_ => True) a
 
         js : Fin m -> List (Fin m)
         js n = filter (<= n) is
