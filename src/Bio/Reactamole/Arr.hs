-- |
-- Module      :  Bio.Reactamole.Arr
-- Copyright   :  (c) DigMP Research Group 2021
-- License     :  MIT
--
-- Arrow definitions for CRNs in Reactamole.
--
-- Sources:
--
--  * /Generalising Monads to Arrows/, by John Hughes,
--    /Science of Computer Programming/ 37, pp67-111, May 2000.
--
--  * /A New Notation for Arrows/, by Ross Paterson, in /ICFP 2001/,
--    Firenze, Italy, pp229-240.
--
--  * [Control.Arrow package.](https://hackage.haskell.org/package/base-4.15.0.0/docs/Control-Arrow.html)

module Bio.Reactamole.Arr
  ( -- * Arrows for CRNs
    (<<<)
  , (>>>)
  , (***)
  , (&&&)
  , first
  , second
  , loop
  --, loop'

    -- * Species Instantiation
  , instCRN

  --, derivativeOfTerm
  --, derivativeOfEq
  --, dropSg
  ) where

import Bio.Reactamole.Core

--------------------------------------------------------------------------------

-- | Right-to-left composition.
(<<<) :: CRN b c -> CRN a b -> CRN a c
(CRN f) <<< (CRN g) = CRN (f . g)

-- | Left-to-right composition.
(>>>) :: CRN a b -> CRN b c -> CRN a c
f >>> g = g <<< f

-- | Split the input between the two argument arrows and combine their output.
--
-- Note that this is in general not a functor.
(***) :: CRN a b -> CRN a' b' -> CRN (a, a') (b, b')
(CRN f) *** (CRN g) = CRN $ \s@(Sg _ sys ic) ->
  let n = length sys
      (Sg spf sysf icf) = (dropSg n . f . runCRN proj1) s
      k = length sysf
      (Sg spg sysg icg) = (shiftSg n k . dropSg n . g . runCRN proj2) s
    in Sg (PairS spf spg) (sys ++ sysf ++ sysg) (ic ++ icf ++ icg)

-- | Fanout: send the input to both argument arrows and combine their output.
(&&&) :: CRN a b -> CRN a c -> CRN a (b, c)
f &&& g = dupCRN >>> (f *** g)

-- | Send the first component of the input through the argument arrow, and copy
-- the rest unchanged to the output.
first :: CRN a b -> CRN (a, c) (b, c)
first f = f *** idCRN

-- | A mirror image of 'first'.
second :: CRN a b -> CRN (c, a) (c, b)
second f = idCRN *** f

--------------------------------------------------------------------------------

-- | Loop the second parameter, c, back into itself.
--
-- The c in the resulting CRN is unused and unmodified by the new CRN
-- (essentially the identity).
loop' :: CRN (a,c) (b,c) -> CRN (a,c) (b,c)
loop' f = CRN $ \s@(Sg (PairS _ y) _ _) ->
  let s'@(Sg (PairS _ y') _ _) = runCRN (second instCRN) s
      s''@(Sg (PairS _ y'') _ _) =  runCRN (f >>> second instCRN) s'
      Sg (PairS x _) sys ic = reWireVars s'' $ zip (varsSp y'') (varsSp y')
  in Sg (PairS x y) sys ic

-- | Feed the second component back in on itself.
loop :: HasDefault c => CRN (a, c) (b, c) -> CRN a b
loop f = CRN $ \(Sg x sys ic) ->
  let CRN g = loop' f
      Sg (PairS y _) sys' ic' = g (Sg (PairS x getDefault) sys ic)
  in Sg y sys' ic'

--------------------------------------------------------------------------------

-- | Instantiate all /uninstantiated/ 'Species'.
instCRN :: CRN a a
instCRN = CRN $ \s@(Sg sp _ _) -> case sp of
  NullS     -> s
  TrueS     -> runCRN instBl s
  FalseS    -> runCRN instBl s
  BoolS _ _ -> runCRN instBl s
  RealS _   -> runCRN instRl s
  PairS _ _ -> runCRN (instCRN *** instCRN) s
  Tup3S {}  -> runCRN instTup3 s
  Tup4S {}  -> runCRN instTup4 s
  Tup5S {}  -> runCRN instTup5 s
  CondS {}  -> runCRN instCond s

-- | Instantiate 'BoolS' species.
instBl :: CRN Bool Bool
instBl = CRN $ \(Sg sp sys ic) ->
  let n = length sys
      constEq = [Term 0.0 []]
  in case sp of
    TrueS     -> Sg (BoolS n (n+1)) (sys ++ [constEq, constEq]) (ic ++ [1,0])
    FalseS    -> Sg (BoolS n (n+1)) (sys ++ [constEq, constEq]) (ic ++ [0,1])
    BoolS _ _ -> Sg sp sys ic

-- | Return derivative of the 'Term'.
derivativeOfTerm :: [Equation] -> Term -> Equation
derivativeOfTerm _ (Term _ []) = []
derivativeOfTerm sys (Term c (x:xs)) =
  let x'  = sys !! x
      xs' = derivativeOfTerm sys (Term 1 xs)
  in addEq (multEq x' [Term c xs]) (multEq [Term c [x]] xs')

-- | Return derivative of the 'Equation'.
derivativeOfEq :: [Equation] -> Equation -> Equation
derivativeOfEq sys = concatMap (derivativeOfTerm sys)

-- | Instantiate 'RealS' species.
instRl :: CRN Double Double
instRl = CRN $ \(Sg (RealS eq) sys ic) ->
  case eq of
    [Term 1 [_], Term (-1) [_]] -> Sg (RealS eq) sys ic
    [Term (-1) [_], Term 1 [_]] -> Sg (RealS (reverse eq)) sys ic
    _ -> let x     = length sys
             x'    = x + 1
             eq'   = [Term 1 [x], Term (-1) [x']]
             eqdt  = derivativeOfEq sys eq
             dxdt  = filter (\(Term c _) -> c > 0) eqdt
             dx'dt = filter (\(Term c _) -> c < 0) eqdt
             sys'  = if null dxdt && null dx'dt
                      then sys ++ [dxdt, dx'dt]
                      else sys ++ [dxdt  ++ [Term (-1) [x, x']],
                                   dx'dt ++ [Term (-1) [x, x']]]
             res   = evalEq eq ic
             x0    = if res > 0 then res  else 0
             x'0   = if res < 0 then -res else 0
             ic'   = ic ++ [x0, x'0]
          in Sg (RealS eq') sys' ic'

-- | Instantiate 'Tup3S' species.
instTup3 :: CRN (a,b,c) (a,b,c)
instTup3 = arrSp (\(Tup3S x y z) -> PairS x (PairS y z))
             >>> (instCRN *** instCRN)
             >>> arrSp (\(PairS x (PairS y z)) -> Tup3S x y z)

-- | Instantiate 'Tup4S' species.
instTup4 :: CRN (a,b,c,d) (a,b,c,d)
instTup4 = arrSp (\(Tup4S x y z w) -> PairS x (Tup3S y z w))
             >>> (instCRN *** instCRN)
             >>> arrSp (\(PairS x (Tup3S y z w)) -> Tup4S x y z w)

-- | Instantiate 'Tup5S' species.
instTup5 :: CRN (a,b,c,d,e) (a,b,c,d,e)
instTup5 = arrSp (\(Tup5S x y z w u) -> PairS x (Tup4S y z w u))
             >>> (instCRN *** instCRN)
             >>> arrSp (\(PairS x (Tup4S y z w u)) -> Tup5S x y z w u)

-- | Instantiate 'CondS' Species.
instCond :: CRN (Either a b) (Either a b)
instCond = arrSp (\(CondS x y z) -> PairS x (PairS y z))
             >>> (instCRN *** instCRN)
             >>> arrSp (\(PairS x (PairS y z)) -> CondS x y z)

--------------------------------------------------------------------------------

-- | Drop the first n variables of a signal.
--
-- Remove the first n _equations_ and the first n _initial conditions_ from the
-- initial value problem
--
-- __WARNING__: This does not re-index the variables! Use with caution!
--              Consider using with 'shiftSg' to re-index them.
dropSg :: Int -> Signal a -> Signal a
dropSg n (Sg sp sys ic) = Sg sp (drop n sys) (drop n ic)
