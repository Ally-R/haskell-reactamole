-- |
-- Module      :  Bio.Reactamole.Bool
-- Copyright   :  (c) DigMP Research Group 2021
-- License     :  MIT
--
-- Tools for Bools: Functions for lifting and manipulating Bools in
-- Reactamole.
--
-- Source: NAND gate construction from /Robust Chemical Circuits/, by Ellis,
-- Klinge and Lathrop, /Biosystems/, 2019.

module Bio.Reactamole.Bool
  ( -- * Signal Functions
    -- ** Constant Boolean
    constBl

    -- ** Logic Gates
  , notCRN
  , nandCRN
  , andCRN
  , orCRN
  , norCRN
  , xorCRN
  , xnorCRN

    -- ** Lift CRNs
  , arr1Bl
  , arr2Bl
  , arr3Bl
  , arr4Bl
  , arr5Bl
  --, arrNOT
  --, arrAND
  --, arrOR
  ) where

import Bio.Reactamole.Arr
import Bio.Reactamole.Core

-- | CRN that produces True (/'instantiate'd/) constantly.
constBl :: Bool -> CRN a Bool
constBl b = constCRN (Sg b' [] [])
  where b' = if b then TrueS else FalseS

--------------------------------------------------------------------------------

-- | CRN for logical NOT.
notCRN :: CRN Bool Bool
notCRN = arrSp $ \sp -> case sp of
  TrueS     -> FalseS
  FalseS    -> TrueS
  BoolS x y -> BoolS y x

--------------------------------------------------------------------------------

-- | CRN for logical NAND.
nandCRN :: CRN (Bool, Bool) Bool
nandCRN = CRN $ \(Sg sp sys ic) -> case sp of
  PairS TrueS  TrueS  -> Sg FalseS [] []
  PairS FalseS _      -> Sg TrueS [] []
  PairS _      FalseS -> Sg TrueS [] []
  PairS (BoolS x1 x1') TrueS -> Sg (BoolS x1' x1) sys ic
  PairS TrueS (BoolS x2 x2') -> Sg (BoolS x2' x2) sys ic
  PairS (BoolS x1 x1') (BoolS x2 x2') ->
    let k  = 0.01
        y  = length sys
        y' = y + 1
        rns = [ Rn [y,   y,  y'] [y,   y,  y ] (3*k)
              , Rn [y',  y', y ] [y',  y', y'] (3*k)
              , Rn [x1,  x2, y ] [x1,  x2, y'] k
              , Rn [x1', y'    ] [x1', y     ] k
              , Rn [x2', y'    ] [x2', y     ] k     ]
    in Sg (BoolS y y') (applyRns rns (sys ++ [[], []])) (ic ++ [0, 1])

--------------------------------------------------------------------------------

-- | CRN for logical AND.
andCRN :: CRN (Bool, Bool) Bool
andCRN = nandCRN >>> notCRN

-- | CRN for logical OR.
orCRN :: CRN (Bool, Bool) Bool
orCRN = (notCRN *** notCRN) >>> nandCRN

-- | CRN for logical NOR.
norCRN :: CRN (Bool, Bool) Bool
norCRN = orCRN >>> notCRN

-- | CRN for logical XOR.
--
-- Note that this uses three nandCRNs under the hood. All other gates except
-- xnorCRN use only one.
xorCRN :: CRN (Bool, Bool) Bool
xorCRN = (nandCRN &&& orCRN) >>> andCRN

-- | CRN for logical XNOR.
--
-- Note that this uses three nandCRNs under the hood. All other gates except
-- xorCRN use only one.
xnorCRN :: CRN (Bool, Bool) Bool
xnorCRN = xorCRN >>> notCRN

--------------------------------------------------------------------------------

-- | Produce either an identical or inverse signal depending on the Boolean.
--
-- Lift helper function.
arrNOT :: Bool -> CRN Bool Bool
arrNOT x = if x then idCRN else notCRN

-- | Compile a series of AND and NOT gates into an CRN, to be fed into OR gates
--  that will simulate a Boolean function.
--
-- Lift helper function.
arrAND :: [Bool] -> Species a -> CRN a Bool
arrAND [] _ = constBl False
arrAND [x] (BoolS _ _) = arrNOT x
arrAND [_] _ = error "Species must be a Bool or Tuple of Bools, and no longer than the Bool list"
arrAND (x:xs) (PairS (BoolS _ _) b@(BoolS _ _)) =
  andCRN <<< (arrNOT x *** arrAND xs b)
arrAND (x:xs) (Tup3S (BoolS _ _) b@(BoolS _ _) c@(BoolS _ _)) =
  andCRN <<< (arrNOT x *** arrAND xs (PairS b c)) <<< tup3ToPairCRN
arrAND (x:xs) (Tup4S (BoolS _ _) b@(BoolS _ _) c@(BoolS _ _) d@(BoolS _ _)) =
  andCRN <<< (arrNOT x *** arrAND xs (Tup3S b c d)) <<< tup4ToPairCRN
arrAND (x:xs) (Tup5S (BoolS _ _)
               b@(BoolS _ _) c@(BoolS _ _) d@(BoolS _ _) e@(BoolS _ _)) =
  andCRN <<< (arrNOT x *** arrAND xs (Tup4S b c d e)) <<< tup5ToPairCRN
arrAND (_:_) _ = error "Bool list cannot be longer than Species; Species must be BoolSpecies or a tuple of BoolSpecies"

-- | Compile a series of OR, AND, and NOT gates into an CRN to simulate a Boolean
--  function.
--
-- Lift helper function.
arrOR :: [[Bool]] -> Species a -> CRN a Bool
arrOR [] s = arrAND [] s
arrOR [x] s = arrAND x s
arrOR (x:xs) s = orCRN <<< (arrAND x s *** arrOR xs s) <<< dupCRN

--------------------------------------------------------------------------------

-- | Lift a 1-in 1-out Boolean Haskell function into an 'CRN'.
arr1Bl :: (Bool -> Bool) -> CRN Bool Bool
arr1Bl f = CRN $ \s@(Sg sp _ _) ->
     let bl = [True, False]
         matrix = [(f x, [x]) | x <- bl]
         mFilt  = map snd . filter ((==True) . fst)
         CRN sg = arrOR (mFilt matrix) sp
     in sg s

-- | Lift a 2-in 1-out Boolean Haskell function into an 'CRN'.
arr2Bl :: (Bool -> Bool -> Bool) -> CRN (Bool, Bool) Bool
arr2Bl f = CRN $ \s@(Sg sp _ _) ->
  let bl = [True, False]
      matrix = [(f x y, [x, y]) | x <- bl, y <- bl]
      mFilt  = map snd . filter ((==True) . fst)
      CRN sg = arrOR (mFilt matrix) sp
  in sg s

-- | Lift a 3-in 1-out Boolean Haskell function into an 'CRN'.
arr3Bl :: (Bool -> Bool -> Bool -> Bool) -> CRN (Bool, Bool, Bool) Bool
arr3Bl f = CRN $ \s@(Sg sp _ _) ->
  let bl = [True, False]
      matrix = [(f x y z, [x, y, z]) | x <- bl, y <- bl, z <- bl]
      mFilt  = map snd . filter ((==True) . fst)
      CRN sg = arrOR (mFilt matrix) sp
  in sg s
      
-- | Lift a 4-in 1-out Boolean Haskell function into an 'CRN'.
arr4Bl :: (Bool -> Bool -> Bool -> Bool -> Bool) ->
          CRN (Bool, Bool, Bool, Bool) Bool
arr4Bl f = CRN $ \s@(Sg sp _ _) ->
  let bl = [True, False]
      matrix = [(f x y z w, [x, y, z, w]) |
                x <- bl, y <- bl, z <- bl, w <- bl]
      mFilt  = map snd . filter ((==True) . fst)
      CRN sg = arrOR (mFilt matrix) sp
  in sg s
  
-- | Lift a 5-in 1-out Boolean Haskell function into an 'CRN'.
arr5Bl :: (Bool -> Bool -> Bool -> Bool -> Bool -> Bool) ->
          CRN (Bool, Bool, Bool, Bool, Bool) Bool
arr5Bl f = CRN $ \s@(Sg sp _ _) ->
  let bl = [True, False]
      matrix = [(f x y z w v, [x, y, z, w, v]) |
                x <- bl, y <- bl, z <- bl, w <- bl, v <- bl]
      mFilt  = map snd . filter ((==True) . fst)
      CRN sg = arrOR (mFilt matrix) sp
  in sg s
