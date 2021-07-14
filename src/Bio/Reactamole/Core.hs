{-# LANGUAGE ScopedTypeVariables #-}

-- |
-- Module      :  Bio.Reactamole.Core
-- Copyright   :  (c) TBD
-- License     :  TBD
--
-- Maintainer  :  TBD
-- Stability   :  TBD
-- Portability :  TBD
--
-- Core Reactamole language structure based on a representation of chemical
-- reaction networks (CRNs) as ordinary differential equations (ODEs).

module Bio.Reactamole.Core
  (  -- * Signal Functions
    SF (..)
  , reduceSF

    -- ** Basic SFs
  , idSF
  , fstSF --__TO_DO:__ Change to proj1
  , sndSF --__TO_DO:__ Change to proj2
  , constSF
  , dupSF

    -- ** Tuple/Pair conversions
  , tup3ToPairSF
  , tup4ToPairSF
  , tup5ToPairSF

    -- * Species & Signals
    -- ** Construction
  , Species (..)
  , Signal (..)
  , HasDefault (..)
  , arrSp

    -- ** Basic Signals
  , nullSg
  , pairSg

    -- ** Manipulation
  , varsSp
  , varsSg
  , reWireSg
  --, shiftSp
  , shiftSg
  --, applyFtoSp
  
    -- * Ordinary Differential Equations
    -- ** Construction
  , Variable
  , Term (..)
  , Equation

    -- ** Manipulation
  , showEq
  , addEq
  , multEq
  , evalEq
  , varsEq
  , normalizeEq
  --, reWireVar
  , reWireVars
  --, dropVar
  --, dropVars
  --, dropIndex
  , varTally
  --, shiftVar
  --, shiftTerm
  --, shiftEq
  --, shiftSys

    -- * Reactions
  , Reaction (..)
  , applyRns
  ) where

import Data.List (intercalate, nub, sort)

--------------------------------------------------------------------------------

-- | Variable of an ODE.
--
-- Acts as an index representing a unique variable name.
type Variable = Int

-- | Term of an ODE.
--
-- For example, @Term 5 [1, 1, 2, 1]@ represents \( \frac{dx}{dt} = 5x^3y \), where
-- \( x \) is represented by 'Variable' @1@ and \( y \) by 'Variable' @2@.
data Term = Term { coef :: Double -- ^Coefficient
                 , vars :: [Variable] -- ^Variables
                 }
            deriving (Eq)

instance Show Term where
  show (Term c vs) = sgn ++ term
    where sgn = if c < 0 then "-" else "+"
          term = intercalate "*x" (show (abs c) : map showVar (varTally vs))
          showVar (a, b) =
            if b == 1 then show a else show a ++ "^" ++ show b

instance Ord Term where
  compare (Term _ v1) (Term _ v2) = compare (sort v1) (sort v2)

-- | Tally a list of variables and remove duplicates.
varTally :: [Variable] -> [(Variable, Int)]
varTally vs = [(x, sum [1 | y <- vs, y == x]) | x <- nub vs]

--------------------------------------------------------------------------------

-- | Ordinary differential equation (ODE).
--
-- For example, @Equation [('Term' 3 [1, 2]), ('Term' -1 [1, 1, 3])]@ represents
-- /dx\/dt = +3*x*y -x*x*z/, where /x/ is represented by 'Vairable' @1@, /y/ by
-- 'Variable' @2@, and /z/ by 'Variable' @3@.
type Equation = [Term]

-- | Turn a polynomial equation into a human readable form.
showEq :: Equation -> String
showEq ts = unwords $ map show ts

-- | Add two polynomial equations together.
addEq :: Equation -> Equation -> Equation
addEq = (++)

-- | Multiply two polynomial equations together.
multEq :: Equation -> Equation -> Equation
multEq eq1 eq2 = [multTerm t1 t2 | t1 <- eq1, t2 <- eq2]
  where multTerm (Term c1 v1) (Term c2 v2) = Term (c1*c2) (v1 ++ v2)

-- | Evaluate an equation with initial conditions.
evalEq  :: Equation -> [Double] -> Double
evalEq eq ic = kernel eq 0
  where kernel [] total = total
        kernel ((Term c vs):ts) total =
          let val = c * product [ic !! v | v <- vs]
          in kernel ts (total + val)

--------------------------------------------------------------------------------

-- | Species definition.
data Species a where
  NullS  :: Species () -- ^ Empty species
  TrueS  :: Species Bool -- ^ Boolean True species /(uninstantiated)/
  FalseS :: Species Bool -- ^ Boolean False species /(uninstantiated)/
  BoolS  :: Variable  -> Variable  -> Species Bool
  RealS  :: Equation  -> Species Double
  PairS  :: Species a -> Species b -> Species (a, b)
  Tup3S  :: Species a -> Species b -> Species c -> Species (a, b, c)
  Tup4S  :: Species a -> Species b -> Species c -> Species d ->
    Species (a, b, c, d)
  Tup5S  :: Species a -> Species b -> Species c -> Species d -> Species e ->
    Species (a, b, c, d, e)
  -- | Conditional species
  CondS  :: Species Bool -> Species a -> Species b -> Species (Either a b)

--------------------------------------------------------------------------------

-- | Default value for a 'Species'.
class HasDefault a where
  getDefault :: Species a

instance HasDefault () where
  getDefault = NullS

instance HasDefault Bool where
  getDefault = FalseS

instance HasDefault Double where
  getDefault = RealS []

instance (HasDefault a, HasDefault b) => HasDefault (a,b) where
  getDefault = PairS (getDefault :: Species a) (getDefault :: Species b)

instance (HasDefault a, HasDefault b, HasDefault c)
  => HasDefault (a,b,c) where
  getDefault = Tup3S (getDefault :: Species a)
                     (getDefault :: Species b)
                     (getDefault :: Species c)

instance (HasDefault a, HasDefault b, HasDefault c, HasDefault d)
  => HasDefault (a,b,c,d) where
  getDefault = Tup4S (getDefault :: Species a)
                     (getDefault :: Species b)
                     (getDefault :: Species c)
                     (getDefault :: Species d)

instance (HasDefault a, HasDefault b, HasDefault c, HasDefault d, HasDefault e)
  => HasDefault (a,b,c,d,e) where
  getDefault = Tup5S (getDefault :: Species a)
                     (getDefault :: Species b)
                     (getDefault :: Species c)
                     (getDefault :: Species d)
                     (getDefault :: Species e)

instance (HasDefault a, HasDefault b) => HasDefault (Either a b) where
  getDefault = CondS (getDefault :: Species Bool)
                     (getDefault :: Species a)
                     (getDefault :: Species b)

--------------------------------------------------------------------------------

instance Show (Species a) where
  show sp = case sp of
    NullS           -> "Null"
    TrueS           -> "True"
    FalseS          -> "False"
    BoolS x y       -> concat ["Bool(x", show x, ", x", show y, ")"]
    RealS eq        -> concat ["Real(", showEq eq, ")"]
    CondS x y z     -> concat ["Cond(", show x, show y, show z, ")"]
    PairS x y       -> tupToStr [show x, show y]
    Tup3S x y z     -> tupToStr [show x, show y, show z]
    Tup4S w x y z   -> tupToStr [show w, show x, show y, show z]
    Tup5S v w x y z -> tupToStr [show v, show w, show x, show y, show z]
    where tupToStr xs = concat ["Tuple(", intercalate ", " xs, ")"]

--------------------------------------------------------------------------------

-- | Get all 'Variable's in a 'Species'.
varsSp :: Species a -> [Variable]
varsSp sp = case sp of
  NullS           -> []
  TrueS           -> []
  FalseS          -> []
  BoolS x y       -> [x, y]
  RealS eq        -> varsEq eq
  CondS x y z     -> nub $ concat [varsSp x, varsSp y, varsSp z]
  PairS x y       -> nub $ varsSp x ++ varsSp y
  Tup3S x y z     -> nub $ concat [varsSp x, varsSp y, varsSp z]
  Tup4S w x y z   -> nub $ concat [varsSp w, varsSp x, varsSp y, varsSp z]
  Tup5S v w x y z -> nub $ concat [varsSp v, varsSp w, varsSp x, varsSp y, varsSp z]

-- | Get all 'Variable's in an 'Equation'.
varsEq :: Equation -> [Variable]
varsEq eq = nub $ concatMap (\(Term _ xs) -> xs) eq

-- | Get all 'Variable's in a 'Signal'.
varsSg :: Signal a -> [Variable]
varsSg = varsSp . species

--------------------------------------------------------------------------------

-- | Apply the given function to each variable in the 'Species'.
applyFtoSp :: (Variable -> Variable) -> Species a -> Species a
applyFtoSp f sp = case sp of
  NullS     -> NullS
  TrueS     -> TrueS
  FalseS    -> FalseS
  BoolS x y -> BoolS (f x) (f y)
  RealS eq  -> RealS (eqF eq f)
  PairS x y       -> PairS (g x) (g y)
  Tup3S x y z     -> Tup3S (g x) (g y) (g z)
  Tup4S x y z w   -> Tup4S (g x) (g y) (g z) (g w)
  Tup5S x y z w u -> Tup5S (g x) (g y) (g z) (g w) (g u)
  CondS x y z -> CondS (g x) (g y) (g z)
  where g :: Species a -> Species a
        g s = applyFtoSp f s
        eqF eq fun = [Term c [fun x | x <- v] | (Term c v) <- eq]

--------------------------------------------------------------------------------

-- | A \"typed\" CRN represented by an ODE.
data Signal a = Sg { species :: Species a  -- ^ /Type/ of the CRN
                   , system  :: [Equation] -- ^ System of equations
                   , init    :: [Double] -- ^ Initial conditions
                   }

-- | A signal function, or CRN transformer.
newtype SF a b = SF { runSF :: Signal a -> Signal b }

instance Show (Signal a) where
  show (Sg sp sys ic) = intercalate "\n" (reverse (map fst eqns)) ++ "\n"
                        ++ intercalate ", " (reverse (map snd eqns))
                        ++ "\nSpecies: " ++ show sp
    where eqns = foldl (\eq (a, b) ->
                          ("dx" ++ show (length eq) ++ "/dt = " ++ showEq a,
                            "x"  ++ show (length eq) ++ "(0) = " ++ show b)
                          : eq)
                        []
                        (zip sys ic)

--------------------------------------------------------------------------------

-- | Identity: emit the input signal unaltered.
idSF :: SF a a
idSF = SF id

-- | The null signal.
nullSg :: Signal ()
nullSg = Sg NullS [] []

-- | Create a Pair signal.
pairSg :: Signal a -> Signal b -> Signal (a, b)
pairSg (Sg sp1 sys1 ic1) (Sg sp2 sys2 ic2) = Sg sp3 sys3 ic3
  where
    n    = length sys1
    sp3  = PairS sp1 (shiftSp 0 n sp2)
    sys3 = sys1 ++ map (shiftEq 0 n) sys2
    ic3  = ic1 ++ ic2

-- | Lift a pure species function to an SF.
--
-- You can regard the species function as a "reindexing" the structure of the
-- species object without actually modifying the underlying system of equations
-- or initial conditions. The resulting SF is really an "adapter" that makes two
-- other SFs compatible.
arrSp :: (Species a -> Species b) -> SF a b
arrSp f = SF $ \(Sg sp sys ic) -> Sg (f sp) sys ic

-- | Extract the first element of a signal.
fstSF :: SF (a, b) a --__TO_DO:__ Change to proj1
fstSF = arrSp $ \(PairS x _) -> x

-- | Extract the second element of a signal.
sndSF :: SF (a, b) b --__TO_DO:__ Change to proj2
sndSF = arrSp $ \(PairS _ y) -> y

--------------------------------------------------------------------------------

-- | Create an SF that ignores its input, producing the constant signal
-- s regardless.
constSF :: Signal a -> SF b a
constSF s = SF (const s)

-- | Duplicate: emit a pair of the input signal.
dupSF :: SF a (a, a)
dupSF = arrSp $ \x -> PairS x x

--------------------------------------------------------------------------------

-- | Convert a 'Tup3S' signal into a 'PairS' signal.
tup3ToPairSF :: SF (a, b, c) (a, (b, c))
tup3ToPairSF = arrSp $ \(Tup3S x y z) -> PairS x (PairS y z)

-- | Convert a 'Tup4S' signal into a 'PairS' signal.
tup4ToPairSF :: SF (a, b, c, d) (a, (b, c, d))
tup4ToPairSF = arrSp $ \(Tup4S x y z w) -> PairS x (Tup3S y z w)

-- | Convert a 'Tup5S' signal into a 'PairS' signal.
tup5ToPairSF :: SF (a, b, c, d, e) (a, (b, c, d, e))
tup5ToPairSF = arrSp $ \(Tup5S x y z w u) -> PairS x (Tup4S y z w u)

--------------------------------------------------------------------------------

-- | Rewire second 'Signal' to first 'Signal', looping the second input
-- 'Species' back into itself.
reWireSg :: Signal (a, c) -> Signal (b, c) -> Signal (b, c)
reWireSg (Sg (PairS _ c) _ _) s@(Sg (PairS _ d) _ _) =
  reWireVars s $ zip (varsSp c) (varsSp d)

-- | Rewire multiple variables.
reWireVars :: Signal a -> [(Variable, Variable)] -> Signal a
reWireVars s = foldr reWireVar s . modTups
  where modTup (x,y) (a,b) =
          if a == b then (x,y) else (shift x a b, shift y a b)
        modTups tups =
          foldr (\x acc -> (foldr (\ac' x' -> modTup x' ac') x acc):acc) [] tups
        shift x a b = if x == b
                  then (if a > b then a-1 else a)
                  else (if x > b then x-1 else x)

-- | Rewire 'Variable' d to 'Variable' c.
reWireVar :: (Variable, Variable) ->  Signal a -> Signal a
reWireVar (c, d) sg@(Sg sp sys ic) =
  if c == d then sg else dropVar (Sg (f sp) sys' ic) d
  where sys' = [[Term con [replaceVar c d v | v <- var] | (Term con var) <- eq]
                | eq <- sys]
        replaceVar :: Variable -> Variable -> Variable -> Variable
        replaceVar cv dv x = if x == dv then cv else x
        f = applyFtoSp (replaceVar c d)

--------------------------------------------------------------------------------

-- | A reaction in a CRN.
data Reaction = Rn { reactants :: [Variable] -- ^ Reactants
                   , products  :: [Variable] -- ^ Products
                   , constant  :: Double -- ^ Rate constant
                   }
                deriving Eq

-- | Reaction displayed similarly to: @x1 + x2 --{5}-> x3 + x3@
instance Show Reaction where
  show (Rn r p k) = toSt r ++ " --{" ++ show k ++ "}-> " ++ toSt p
    where toSt = intercalate " + " . map (("x" ++) . show)

instance Ord Reaction where
  compare (Rn r1 _ k1) (Rn r2 _ k2) = compare (sort r1,k1) (sort r2,k2)

-- | Update a system of 'Equation's to include the effects of a list of
-- 'Reaction's.
applyRns :: [Reaction] -> [Equation] -> [Equation]
applyRns rns sys = map (normalizeEq . extendEq) [0..length sys - 1]
  where extendEq n = sys !! n ++ map (toTerm n) rns
        toTerm x (Rn r p k) = Term (k * fromIntegral (count x p - count x r)) r
        count x = length . filter (==x)

-- | Simplify the 'Equation' so that all terms are combined and sorted and all
-- zero terms are dropped.
normalizeEq :: Equation -> Equation
normalizeEq = filter notZero . combineTerms . sort . map sortVars
  where notZero (Term c _) = c /= 0
        sortVars (Term c vs) = Term c (sort vs)
        combineTerms [] = []
        combineTerms [t] = [t]
        combineTerms ( t1@(Term c1 v1) : t2@(Term c2 v2) : ts)
          | v1 == v2  = combineTerms $ Term (c1 + c2) v1 : ts
          | otherwise = t1 : combineTerms (t2 : ts)

--------------------------------------------------------------------------------

-- | Remove unnecessary variables in the signal.
--
-- This is accomplished by finding the smallest "closed sub-CRN" in the
-- signal that encapsulates all of the variables in the output.
reduceSF :: SF a a
reduceSF = SF $ \s@(Sg sp sys _) ->
  let n = length sys
      deps = map (\i -> varsEq (sys !! i)) [0..n-1]
      f vs = let vs' = sort . nub $ concat (vs : [ deps !! v | v <- vs])
             in if vs == vs' then vs else f vs'
      used = f (sort $ varsSp sp)
      notUsed = filter (not . (`elem` used)) [0..n-1]
  in dropVars s notUsed

-- | Remove all of the variables from the signal.
--
-- WARNING: This function simply removes the relevant equations and
--          initial conditions from the signal. It is possible that
--          the variables are used in other equations.
dropVars :: Signal a -> [Variable] -> Signal a
dropVars s = foldl dropVar s . sort

-- | Remove the given variable from the signal.
--
-- WARNING: This function simply removes the relevant equations and
--          initial conditions from the signal. It is possible that
--          the variables are used in other equations.
dropVar :: Signal a -> Variable -> Signal a
dropVar s v =
  let Sg sp sys ic = shiftSg v (-1) s
      sys' = dropIndex v sys
      ic'  = dropIndex v ic
  in Sg sp sys' ic'

-- | Remove the element at index n from a list.
dropIndex :: Int -> [a] -> [a]
dropIndex n xs = lft ++ rgt
  where (lft, _:rgt) = splitAt n xs

-- | Shift variables whose index is @>= n@ by @+k@ in the species and system of
-- equations of the signal.
--
-- For example, @shiftSg 4 3 s@  will shift signal @s@ by modifying all
-- variable indices starting at @4..@ by @+3@.
--
-- __WARNING__: This does not verify the number of equations is the same
-- as the indices of the variables! Use with caution!
shiftSg :: Int -> Int -> Signal a -> Signal a
shiftSg n k (Sg sp sys ic) = Sg (shiftSp n k sp) (shiftSys n k sys) ic

-- | Shift variables whose index is @>= n@ by @+k@.
shiftVar :: Int -> Int -> Variable -> Variable
shiftVar n k var = if var < n then var else var + k

-- | Shift variables whose index is @>= n@ by @+k@.
shiftTerm :: Int -> Int -> Term -> Term
shiftTerm n k (Term c v) = Term c (map (shiftVar n k) v)

-- | Shift variables whose index is @>= n@ by @+k@.
shiftEq :: Int -> Int -> Equation -> Equation
shiftEq n k = map (shiftTerm n k)

-- | Shift variables whose index is @>= n@ by @+k@.
shiftSys :: Int -> Int -> [Equation] -> [Equation]
shiftSys n k = map (shiftEq n k)

-- | Shift variables whose index is @>= n@ by @+k@.
shiftSp :: Int -> Int -> Species a -> Species a
shiftSp n k sp = case sp of
  NullS     -> NullS
  TrueS     -> TrueS
  FalseS    -> FalseS
  BoolS x y -> BoolS (shiftVar n k x) (shiftVar n k y)
  RealS eq  -> RealS (shiftEq n k eq)
  PairS x y       -> PairS (f x) (f y)
  Tup3S x y z     -> Tup3S (f x) (f y) (f z)
  Tup4S x y z w   -> Tup4S (f x) (f y) (f z) (f w)
  Tup5S x y z w u -> Tup5S (f x) (f y) (f z) (f w) (f u)
  CondS x y z -> CondS (f x) (f y) (f z)
  where f :: Species a -> Species a
        f = shiftSp n k
