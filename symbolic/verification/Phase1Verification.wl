(* ::Package:: *)

(* ============================================================ *)
(* Phase I: dimensional-to-nondimensional verification           *)
(* ============================================================ *)

(* Purpose:
   Verify the Kelvin--Voigt dispersion relation transformation

     D(gamma) -> DStarLambda(s).

   Mathematical statement tested:

     DStarLambda[] == rhoT/(2 k^2) * DimensionalDispersion[]

   after applying NondimensionalizationRules[] and under Phase1Assumptions[].

   Dependencies:
     None.

   Authoritative objects exported by this file:
     - DimensionalDispersion[]
     - NondimensionalizationRules[]
     - DStarLambda[]
     - ScaledDimensionalD[]
     - RunPhase1Verification[]

   Expected pass/fail behavior:
     - PASS requires every component-wise mapping check and the full
       expression equivalence check to verify exactly.

   This file intentionally stops at D(gamma) -> DStarLambda(s). It does not
   derive the polynomial representation and does not alter the physics.
*)

ClearAll[
  gamma, k, rho1, rho2, mu1, mu2, G1, G2,
  rhoT, muT, s, Arho, Amu, AG, Lambda,
  Phase1Symbols, DimensionalDomainAssumptions,
  NondimensionalDomainAssumptions, Phase1Assumptions,
  DimensionalEta1, DimensionalEta2,
  DimensionalDenominators, DimensionalRadicands, DimensionalRadicals,
  DimensionalDispersion,
  OriginalExpressionDomainConditions,
  NondimensionalizationRules,
  Eta1Star, Eta2Star, NondimensionalRadicands, Q1Star, Q2Star,
  NondimensionalDenominators,
  DStarLambda,
  ScalingPrefactor, ScaledDimensionalD,
  EtaScalingDifferences, RadicandMappingDifferences,
  RadicalMappingDifferences, DenominatorScalingDifferences,
  ConstantTermDifference, ReciprocalTermDifference,
  Phase1Difference, VerifyPhase1Equivalence,
  Phase1VerificationTable,
  RunPhase1Verification
];

Phase1Symbols[] := {
  gamma, k, rho1, rho2, mu1, mu2, G1, G2,
  rhoT, muT, s, Arho, Amu, AG, Lambda
};

DimensionalDomainAssumptions[] := And[
  Element[gamma, Complexes],
  Element[{k, rho1, rho2, mu1, mu2, G1, G2}, Reals],
  gamma != 0,
  k > 0,
  rho1 > 0,
  rho2 > 0,
  mu1 > 0,
  mu2 > 0,
  G1 >= 0,
  G2 >= 0,
  G1 + G2 > 0
];

NondimensionalDomainAssumptions[] := And[
  Element[{k, rhoT, muT, Arho, Amu, AG, Lambda}, Reals],
  Element[s, Complexes],
  k > 0,
  rhoT > 0,
  muT > 0,
  Lambda > 0,
  -1 < Arho < 1,
  -1 < Amu < 1,
  -1 <= AG <= 1,
  s != 0
];

Phase1Assumptions[] := NondimensionalDomainAssumptions[];

DimensionalEta1[] := mu1 + G1/gamma;

DimensionalEta2[] := mu2 + G2/gamma;

DimensionalDenominators[] := {
  DimensionalEta1[] + DimensionalEta2[]*
    Sqrt[1 + rho2*gamma/(DimensionalEta2[]*k^2)],
  DimensionalEta2[] + DimensionalEta1[]*
    Sqrt[1 + rho1*gamma/(DimensionalEta1[]*k^2)]
};

DimensionalRadicands[] := {
  1 + rho2*gamma/(DimensionalEta2[]*k^2),
  1 + rho1*gamma/(DimensionalEta1[]*k^2)
};

DimensionalRadicals[] := {
  Sqrt[DimensionalRadicands[][[1]]],
  Sqrt[DimensionalRadicands[][[2]]]
};

DimensionalDispersion[] :=
  gamma*(1/DimensionalDenominators[][[1]] +
      1/DimensionalDenominators[][[2]]) +
    4*k^2/(rho1 + rho2);

OriginalExpressionDomainConditions[] := And[
  gamma != 0,
  k != 0,
  rho1 + rho2 != 0,
  DimensionalEta1[] != 0,
  DimensionalEta2[] != 0,
  DimensionalDenominators[][[1]] != 0,
  DimensionalDenominators[][[2]] != 0
];

NondimensionalizationRules[] := {
  gamma -> s*muT*k^2/rhoT,
  rho1 -> rhoT*(1 + Arho)/2,
  rho2 -> rhoT*(1 - Arho)/2,
  mu1 -> muT*(1 + Amu)/2,
  mu2 -> muT*(1 - Amu)/2,
  G1 -> Lambda*muT^2*k^2*(1 + AG)/(2*rhoT),
  G2 -> Lambda*muT^2*k^2*(1 - AG)/(2*rhoT)
};

Eta1Star[] := (1 + Amu) + (Lambda/s)*(1 + AG);

Eta2Star[] := (1 - Amu) + (Lambda/s)*(1 - AG);

NondimensionalRadicands[] := {
  1 + ((1 - Arho)*s)/Eta2Star[],
  1 + ((1 + Arho)*s)/Eta1Star[]
};

Q1Star[] := Sqrt[1 + ((1 + Arho)*s)/Eta1Star[]];

Q2Star[] := Sqrt[1 + ((1 - Arho)*s)/Eta2Star[]];

NondimensionalDenominators[] := {
  Eta1Star[] + Eta2Star[]*Q2Star[],
  Eta2Star[] + Eta1Star[]*Q1Star[]
};

DStarLambda[] :=
  s*(1/(Eta1Star[] + Eta2Star[]*Q2Star[]) +
      1/(Eta2Star[] + Eta1Star[]*Q1Star[])) + 2;

ScalingPrefactor[] := rhoT/(2*k^2);

ScaledDimensionalD[] :=
  FullSimplify[
    ScalingPrefactor[]*
      (DimensionalDispersion[] /. NondimensionalizationRules[]),
    Phase1Assumptions[]
  ];

EtaScalingDifferences[] :=
  FullSimplify[
    {
      (DimensionalEta1[] /. NondimensionalizationRules[])/(muT/2) -
        Eta1Star[],
      (DimensionalEta2[] /. NondimensionalizationRules[])/(muT/2) -
        Eta2Star[]
    },
    Phase1Assumptions[]
  ];

RadicandMappingDifferences[] :=
  FullSimplify[
    (DimensionalRadicands[] /. NondimensionalizationRules[]) -
      NondimensionalRadicands[],
    Phase1Assumptions[]
  ];

RadicalMappingDifferences[] :=
  FullSimplify[
    (DimensionalRadicals[] /. NondimensionalizationRules[]) -
      {Q2Star[], Q1Star[]},
    Phase1Assumptions[]
  ];

DenominatorScalingDifferences[] :=
  FullSimplify[
    {
      (DimensionalDenominators[][[1]] /. NondimensionalizationRules[])/(muT/2) -
        NondimensionalDenominators[][[1]],
      (DimensionalDenominators[][[2]] /. NondimensionalizationRules[])/(muT/2) -
        NondimensionalDenominators[][[2]]
    },
    Phase1Assumptions[]
  ];

ConstantTermDifference[] :=
  FullSimplify[
    ScalingPrefactor[]*
      ((4*k^2/(rho1 + rho2)) /. NondimensionalizationRules[]) - 2,
    Phase1Assumptions[]
  ];

ReciprocalTermDifference[] :=
  FullSimplify[
    ScalingPrefactor[]*
      (gamma*(1/DimensionalDenominators[][[1]] +
           1/DimensionalDenominators[][[2]]) /.
        NondimensionalizationRules[]) -
      s*(1/NondimensionalDenominators[][[1]] +
          1/NondimensionalDenominators[][[2]]),
    Phase1Assumptions[]
  ];

Phase1Difference[] :=
  FullSimplify[
    ScaledDimensionalD[] - DStarLambda[],
    Phase1Assumptions[]
  ];

VerifyPhase1Equivalence[] := TrueQ[Phase1Difference[] === 0];

Phase1VerificationTable[] := {
  <|"Check" -> "eta scaling", "Expected" -> "{0, 0}",
    "Result" -> EtaScalingDifferences[],
    "Passed" -> TrueQ[EtaScalingDifferences[] === {0, 0}]|>,
  <|"Check" -> "radicand mapping", "Expected" -> "{0, 0}",
    "Result" -> RadicandMappingDifferences[],
    "Passed" -> TrueQ[RadicandMappingDifferences[] === {0, 0}]|>,
  <|"Check" -> "principal radical mapping", "Expected" -> "{0, 0}",
    "Result" -> RadicalMappingDifferences[],
    "Passed" -> TrueQ[RadicalMappingDifferences[] === {0, 0}]|>,
  <|"Check" -> "denominator scaling", "Expected" -> "{0, 0}",
    "Result" -> DenominatorScalingDifferences[],
    "Passed" -> TrueQ[DenominatorScalingDifferences[] === {0, 0}]|>,
  <|"Check" -> "constant term scaling", "Expected" -> "0",
    "Result" -> ConstantTermDifference[],
    "Passed" -> TrueQ[ConstantTermDifference[] === 0]|>,
  <|"Check" -> "reciprocal term scaling", "Expected" -> "0",
    "Result" -> ReciprocalTermDifference[],
    "Passed" -> TrueQ[ReciprocalTermDifference[] === 0]|>,
  <|"Check" -> "full expression equivalence", "Expected" -> "0",
    "Result" -> Phase1Difference[],
    "Passed" -> VerifyPhase1Equivalence[]|>
};

RunPhase1Verification[] := <|
  "Phase" -> "Phase I: D(gamma) -> DStarLambda(s)",
  "DimensionalDomainAssumptions" -> DimensionalDomainAssumptions[],
  "NondimensionalDomainAssumptions" -> NondimensionalDomainAssumptions[],
  "Assumptions" -> Phase1Assumptions[],
  "OriginalExpressionDomainConditions" -> OriginalExpressionDomainConditions[],
  "ScalingPrefactor" -> ScalingPrefactor[],
  "DimensionalDispersion" -> DimensionalDispersion[],
  "DimensionalDenominators" -> DimensionalDenominators[],
  "DimensionalRadicands" -> DimensionalRadicands[],
  "DimensionalRadicals" -> DimensionalRadicals[],
  "NondimensionalizationRules" -> NondimensionalizationRules[],
  "EtaScalingDifferences" -> EtaScalingDifferences[],
  "RadicandMappingDifferences" -> RadicandMappingDifferences[],
  "RadicalMappingDifferences" -> RadicalMappingDifferences[],
  "DenominatorScalingDifferences" -> DenominatorScalingDifferences[],
  "ConstantTermDifference" -> ConstantTermDifference[],
  "ReciprocalTermDifference" -> ReciprocalTermDifference[],
  "VerificationTable" -> Phase1VerificationTable[],
  "DStarLambda" -> DStarLambda[],
  "ScaledDimensionalD" -> ScaledDimensionalD[],
  "Difference" -> Phase1Difference[],
  "EquivalenceVerified" -> VerifyPhase1Equivalence[],
  "TransformationType" ->
    "Equivalence up to nonzero multiplicative factor rhoT/(2 k^2), assuming k > 0 and rhoT > 0.",
  "NotYetPerformed" -> {
    "Polynomial derivation",
    "Radical elimination",
    "Denominator clearing",
    "Admissible-domain construction",
    "Spurious-root classification"
  }
|>;
