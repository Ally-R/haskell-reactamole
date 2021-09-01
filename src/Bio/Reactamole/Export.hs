-- |
-- Module      :  Bio.Reactamole.Export
-- Copyright   :  (c) DigMP Research Group 2021
-- License     :  MIT
--
-- A collection of useful functions for cleaning up, printing, and converting
-- between ODEs and CRNs.

module Bio.Reactamole.Export where

import Bio.Reactamole.Core
import Bio.Reactamole.Arr
import Data.List (delete, groupBy, sort, nub, intercalate)

data IVP a b  = IVP (Species a) (Species b) [Equation] [Double]
data Rxns a b = Rxns (Species a) (Species b) [Reaction] [Double]

showIn,showOut :: Species a -> String
showIn sp  = "INPUT:\n  "  ++ show sp ++ "\n"
showOut sp = "OUTPUT:\n  " ++ show sp ++ "\n"

showSys :: [Variable] -> [Equation] -> String
showSys inputs sys = "EQUATIONS:\n" ++ concat (zipWith displayEq [0..length sys-1] sys)
  where displayEq x eq =
          if x `elem` inputs
            then ""
            else "  dx" ++ show x ++ "/dt = " ++ show eq ++ "\n"

showICs :: [Variable] -> [Double] -> String
showICs inputs ic = "INITIAL CONDITIONS:\n" ++ concat (zipWith showIC [0..length ic - 1] ic)
  where showIC x x0 =
          if x `elem` inputs
            then ""
            else "  x" ++ show x ++ "(0) = " ++ show x0 ++ "\n"

showRns :: [Reaction] -> String
showRns rns = "REACTIONS:\n" ++ intercalate "\n" (map (("  " ++) . show) rns) ++ "\n"

instance Show (IVP a b) where
  show (IVP spa spb sys ic) =
    intercalate "\n" [showIn spa, showOut spb, showSys inputs sys, showICs inputs ic]
    where inputs = varsSp spa

instance Show (Rxns a b) where
  show (Rxns spa spb rns ic) =
    intercalate "\n" [showIn spa, showOut spb, showRns rns, showICs (varsSp spa) ic]

-- | Extract the initial value problem from the signal function.
toIVP :: HasDefault a => CRN a b -> IVP a b
toIVP f = IVP inSp outSp (map normalizeEq sys) ic
  where s@(Sg inSp _ _) = runCRN instCRN (Sg getDefault [] [])
        f' = f >>> instCRN >>> reduceCRN
        Sg outSp sys ic = runCRN f' s

termToRn :: Variable -> Term -> Reaction
termToRn x (Term c ys) =
  if c > 0 then Rn ys (x:ys) c
           else Rn ys (delete x ys) (-c)

eqsToRns :: [Equation] -> [Reaction]
eqsToRns sys = sort (map (foldr1 combineRns) (groupRns rns))
  where rns = concat $ zipWith (map . termToRn) [0..length sys - 1] sys

groupRns :: [Reaction] -> [[Reaction]]
groupRns = groupBy f . sort
  where f (Rn r1 _ k1) (Rn r2 _ k2) = sort r1 == sort r2 && k1 == k2

combineRns :: Reaction -> Reaction -> Reaction
combineRns (Rn r p1 k) (Rn _ p2 _) = Rn r p k
  where tr  = varTally r
        tr' = map (\(x,n) -> (x,-n)) tr
        tp1 = varTally p1
        tp2 = varTally p2
        ne1 = addTallies tp1 tr'
        ne2 = addTallies tp2 tr'
        tp  = addTallies tr (addTallies ne1 ne2)
        p   = concatMap (\(x,n) -> replicate n x) tp

addTallies :: [(Variable,Int)] -> [(Variable,Int)] -> [(Variable,Int)]
addTallies xs ys = map total (nub (map fst zs))
  where zs = xs ++ ys
        total z = (z,sum (map snd (filter ((==z) . fst) zs)))

toRxns :: HasDefault a => CRN a b -> Rxns a b
toRxns f = Rxns spa spb (eqsToRns sys) ic
  where IVP spa spb sys ic = toIVP f
