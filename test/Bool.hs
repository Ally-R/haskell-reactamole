module Bool where

import Test.Tasty
import Test.Tasty.HUnit

import Bio.Reactamole

truishConc :: [Double] -> Bool
truishConc [_, x] = x > eps
  where
    eps = 0.001
truishConc _      = error "(truishConc) given more than two concentrations"

signalBoolTestCase :: String -> Signal Bool -> Bool -> TestTree
signalBoolTestCase desc sg expected = testCase desc test
  where
    time = 2000.0
    test = do
      table <- calculateSignalTable time sg
      let results = last table
      truishConc results == expected @? "Expected: " ++ show expected ++ ", received: " ++ show results

allTestsUn :: String -> CRN Bool Bool -> [Bool] -> TestTree
allTestsUn desc sf exps =
  testGroup desc $ zipWith3 mkCase allInputs [1..2] exps
    where
      allInputs = [trueSg, falseSg]
      mkCase initial n expected =
        signalBoolTestCase (show (n :: Integer))
                           (runCRN (sf >>> instBl) initial)
                           expected

allTestsUnI :: String -> CRN Bool Bool -> [Bool] -> TestTree
allTestsUnI desc sf exps =
  testGroup desc $ zipWith3 mkCase allInputs [1..2] exps
    where
      allInputs = [trueSg, falseSg]
      mkCase b n expected =
        signalBoolTestCase (show (n :: Integer))
                           (runCRN (constCRN b >>> instBl >>> sf) nullSg)
                           expected

allTestsBin :: String -> CRN (Bool, Bool) Bool -> [Bool] -> TestTree
allTestsBin desc sf exps =
  testGroup desc $ zipWith3 mkCase allInputs [1..4] exps
  where
    allInputs =
      map (uncurry pairSg) [ (trueSg, trueSg)
                           , (trueSg, falseSg)
                           , (falseSg, trueSg)
                           , (falseSg, falseSg)
                           ]
    mkCase initial n expected =
      signalBoolTestCase (show (n :: Integer))
                         (runCRN (sf >>> instBl) initial)
                         expected

allTestsBinI :: String -> CRN (Bool, Bool) Bool -> [Bool] -> TestTree
allTestsBinI desc sf exps =
  testGroup desc $ zipWith3 mkCase allInputs [1..4] exps
  where
    allInputs = [ (trueSg, trueSg)
                , (trueSg, falseSg)
                , (falseSg, trueSg)
                , (falseSg, falseSg)
                ]
    mkCase (b1, b2) n expected =
      signalBoolTestCase (show (n :: Integer))
                         (runCRN ((constCRN b1 >>> instBl) &&& (constCRN b2 >>> instBl) >>> sf) nullSg)
                         expected

--------------------------------------------------------------------------------

boolTests :: TestTree
boolTests = testGroup "Boolean Operators"
  [ unaryEarly
  , binaryEarly
  , unaryLate
  , binaryLate
  ]

unaryEarly :: TestTree
unaryEarly = testGroup "Unary (Early Instantiation)"
  [ allTestsUnI "idCRN"  idCRN  [True, False]
  , allTestsUnI "notCRN" notCRN [False, True]
  ]

binaryEarly :: TestTree
binaryEarly = testGroup "Binary (Early Instantiation)"
  [ allTestsBinI "nandCRN" nandCRN [False, True, True, True]
  , allTestsBinI "andCRN"  andCRN  [True, False, False, False]
  , allTestsBinI "orCRN"   orCRN   [True, True, True, False]
  , allTestsBinI "norCRN"  norCRN  [False, False, False, True]
  , allTestsBinI "xorCRN"  xorCRN  [False, True, True, False]
  ]

unaryLate :: TestTree
unaryLate = testGroup "Unary (Late Instantiation)"
  [ allTestsUn "idCRN"  idCRN  [True, False]
  , allTestsUn "notCRN" notCRN [False, True]
  ]

binaryLate :: TestTree
binaryLate = testGroup "Binary (Late Instantiation)"
  [ allTestsBin "nandCRN" nandCRN [False, True, True, True]
  , allTestsBin "andCRN"  andCRN  [True, False, False, False]
  , allTestsBin "orCRN"   orCRN   [True, True, True, False]
  , allTestsBin "norCRN"  norCRN  [False, False, False, True]
  , allTestsBin "xorCRN"  xorCRN  [False, True, True, False]
  ]
