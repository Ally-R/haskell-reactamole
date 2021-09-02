-- |
-- Module      :  Bio.Reactamole.Matlab
-- Copyright   :  (c) DigMP Research Group 2021
-- License     :  MIT
--
-- Backend to simulate CRNs outputted from Reactamole in MATLAB.

module Bio.Reactamole.Matlab
  ( printSignal
  , printSignalAll
  , displaySignal
  , displaySignalAll
  , saveSignal
  , saveSignalAll
  , calculateSignalTable
  , spStr

  --, sysStr
  --, sysToDecl
  --, simToDecl
  --, spToDecl
  --, allVarsToDecl
  --, waitAndQuit
  --, saveAndQuit
  --, displayTableAndQuit
  --, sgToTable
  --, sgToScript
  --, sgAllToScript
  --, callMatlabWith
  --, readMatlabWith
  )
  where

import Bio.Reactamole.Core

import Data.List (intercalate)
import qualified Data.Text as T
import System.Process (callProcess, readProcess)

--------------------------------------------------------------------------------

sysStr :: [Equation] -> String
sysStr sys = "[" ++ intercalate "; " (map eqStr sys) ++ "]"
  where eqStr [] = "0"
        eqStr eq = intercalate " + " (map termStr eq)
        termStr (Term c vs) = intercalate "*" (show c : map varStr vs)
        varStr v = "x(" ++ show (v + 1) ++ ")"

spStr :: Species a -> [String]
spStr s = case s of
  NullS     -> []
  TrueS     -> ["t, ones(size(t))"]
  FalseS    -> ["t, ones(size(t))*0"]
  BoolS x _ -> ["t, " ++ varStr x]
  RealS eq  -> ["t, " ++ eqStr eq]
  PairS x y -> spStr x ++ spStr y
  Tup3S x y z -> concat [spStr x, spStr y, spStr z]
  Tup4S x y z w -> concat [spStr x, spStr y, spStr z, spStr w]
  Tup5S x y z w u -> concat [spStr x, spStr y, spStr z, spStr w, spStr u]
  CondS x y z -> concat [spStr x, spStr y, spStr z]
  where eqStr [] = "ones(size(t))*0"
        eqStr eq = intercalate " + " (map termStr eq)
        termStr (Term c []) = "ones(size(t))*" ++ show c
        termStr (Term c vs) = intercalate "*" (show c : map varStr vs)
        varStr x = "x(:," ++ show (x + 1) ++ ")"

--------------------------------------------------------------------------------

sysToDecl :: [Equation] -> String
sysToDecl []  = ""
sysToDecl sys = "sys = @(t,x) [" ++ sysStr sys ++ "];"

simToDecl :: Double -> [Double] -> String
simToDecl _ [] = ""
simToDecl t ic = "[t,x] = ode45(sys, [0 " ++ show t ++ "], " ++ initStr
  where initStr = "[" ++ intercalate "; " (map show ic) ++ "]);"

spToDecl :: Species a -> String
spToDecl sp = "plot(" ++ intercalate ", " (spStr sp) ++ ");"

allVarsToDecl :: [Equation] -> String
allVarsToDecl sys =  "plot(" ++ plotStr (length sys) ++ ");"
  where plotStr n = intercalate ", " (map seriesStr [1..n])
        seriesStr x = "t, x(:," ++ show x ++ ")"

--------------------------------------------------------------------------------

type MatlabScript = [String]

waitAndQuit :: MatlabScript
waitAndQuit = ["f = gcf;", "uiwait(f);", "quit;"]

saveAndQuit :: FilePath -> MatlabScript
saveAndQuit path = ["f = gcf;", "hgsave(f,'" ++ path ++ "');", "quit;"]

displayTableAndQuit :: Species a -> MatlabScript
displayTableAndQuit sp = ["fprintf(" ++ fmt ++ "," ++ intercalate "," (spStr sp) ++ ");", "quit;"]
  where
    n   = (length $ spStr sp) + 1   -- N.B., 1 extra column for t
    fmt = "\"" ++ (intercalate " " $ replicate n "%f") ++ "\\n\""

-- | Produce a Matlab script that calculates the numeric solution of the signal.
sgToTable :: Double -> Signal a -> MatlabScript
sgToTable t (Sg _ sys ic) = [sysToDecl sys, simToDecl t ic]

-- | Print a signal to three lines of MATLAB code that can be executed
-- immediately to plot the signal
sgToScript :: Double -> Signal a -> MatlabScript
sgToScript t sg@(Sg sp _ _) = sgToTable t sg ++ [spToDecl sp]

-- | Print all the variables of the signal and ignores the species part
-- of the signal. Can be useful for debugging.
sgAllToScript :: Double -> Signal a -> MatlabScript
sgAllToScript t sg@(Sg _ sys _) = sgToTable t sg ++ [allVarsToDecl sys]

--------------------------------------------------------------------------------

-- | 'printSingal t sg' prints a Matlab script to display 'sg' over 't'
-- timesteps.
printSignal :: Double -> Signal a -> IO ()
printSignal t sg = putStrLn $ unwords $ sgToScript t sg

-- | 'printSingalAll t sg' prints a Matlab script to display all the species of
-- 'sg' over 't' timesteps.
printSignalAll :: Double -> Signal a -> IO ()
printSignalAll t sg = putStrLn $ intercalate "\n" $ sgAllToScript t sg

callMatlabWith :: MatlabScript -> IO ()
callMatlabWith script = callProcess "matlab" ["-nosplash", "-batch", unwords script]

readMatlabWith :: MatlabScript -> IO String
readMatlabWith script = readProcess "matlab" ["-nosplash", "-batch", unwords script] ""

-- | 'displaySingal t sg' displays signal 'sg' over 't' time steps in Matlab.
displaySignal :: Double -> Signal a -> IO ()
displaySignal t sg =
  callMatlabWith $ sgToScript t sg ++ waitAndQuit

-- | 'displaySingalAll t sg' displays all of the species found in signal 'sg'
-- over 't' time steps in Matlab.
displaySignalAll :: Double -> Signal a -> IO ()
displaySignalAll t sg =
  callMatlabWith $  sgAllToScript t sg ++ waitAndQuit

-- | 'saveSingal path t sg' saves a plot of 'sg' over 't' time steps to 'path'.
saveSignal :: FilePath -> Double -> Signal a -> IO ()
saveSignal path t sg =
  callMatlabWith $ sgToScript t sg ++ saveAndQuit path

-- | 'saveSingalAll path t sg' saves a plot of species found in signal 'sg'
-- over 't' time steps to 'path'.
saveSignalAll :: FilePath -> Double -> Signal a -> IO ()
saveSignalAll path t sg =
  callMatlabWith $  sgAllToScript t sg ++ saveAndQuit path

-- | 'calculateSignalTable t sg' queries Matlab to produce the numeric solution
-- for 'sg' over time interval 't'. The solution is returned as a table of
-- points where rows correspond to times and columns correspond to entries for
-- each variable of the signal's equation.
calculateSignalTable :: Double -> Signal a -> IO [[Double]]
calculateSignalTable t sg = do
  -- putStrLn $ intercalate "\n" $ sgToTable t sg ++ displayTableAndQuit (species sg)
  str <- readMatlabWith $ sgToTable t sg ++ displayTableAndQuit (species sg)
  let txt = T.strip (T.pack str)
      ls = T.splitOn (T.pack "\n") txt
  pure $ [ map (\w -> read (T.unpack w) :: Double) l | l <- map T.words ls ]
