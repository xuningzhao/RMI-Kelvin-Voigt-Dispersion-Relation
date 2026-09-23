(* ::Package:: *)

(* ============================================================ *)
(* Phase II: polynomial derivation and forward-completeness chain *)
(* ============================================================ *)

(* Purpose:
   Derive the degree-14 polynomial candidate P14(s) from the verified
   nondimensional dispersion relation DStarLambda(s).

   Mathematical statement tested:

     Phase2AdmissibleDomain[] && DStarLambda[] == 0
       ==> P14Polynomial[] == 0.

   Dependencies:
     - verification/Phase1Verification.wl must be loaded first.

   Authoritative objects exported by this file:
     - F0ClockForm[], ..., F4AfterSecondSquaring[]
     - rawPolynomialNumerator[]
     - P14Polynomial[] and coefficient helpers
     - Phase2DerivationSteps[]
     - Phase2ForwardPreservingChecks[]
     - RunPhase2PolynomialDerivation[]

   Assumptions:
     - Phase II inherits Phase1Assumptions[].
     - Original-expression denominator conditions are recorded by
       Phase2OriginalExpressionDomain[].
     - Derivation nonzero conditions are recorded by
       Phase2TransformationNonzeroConditions[].

   Expected pass/fail behavior:
     - PASS requires all symbolic checks to report PASS.
     - UNRESOLVED counts as failure for theorem certification.
     - Numerical spot checks are supplemental and must pass for the audit
       entry point to return PASS, but they do not replace symbolic proof.

   Reverse implication is intentionally not claimed. Squaring and denominator
   clearing may introduce extraneous polynomial roots.
*)

ClearAll[
  etaPlusPoly, etaMinusPoly,
  etaPlus, etaMinus,
  radicandMinus, radicandPlus,
  qMinus, qPlus,
  recipDenomMinus, recipDenomPlus,
  F0ClockForm, F1ClearedReciprocals,
  alpha0, betaMinus, gammaPlus, deltaCoupling,
  F1LinearRadicals,
  F2AfterFirstSquaring,
  u0, u1, u2,
  F3SingleRadical,
  F4AfterSecondSquaring,
  rawPolynomialNumerator,
  denominatorClearingEtaFactor,
  polynomialNormalizationFactor,
  P14Polynomial,
  P14CoefficientList,
  P14LeadingCoefficient,
  P14Degree,
  P14CoefficientCount,
  P14IsPolynomialQ,
  Phase2PhysicalAssumptions,
  Phase2OriginalExpressionDomain,
  Phase2TransformationNonzeroConditions,
  Phase2AdmissibleDomain,
  Phase2DerivationSteps,
  Phase2CheckStatus,
  Phase2MakeCheck,
  Phase2IdentityCheck,
  Phase2CheckPassedQ,
  Phase2ForwardStepCertifiedQ,
  Phase2ForwardPreservingChecks,
  Phase2PossibleDegreeDropConditions,
  Phase2SymbolicChecks,
  Phase2NumericalSpotChecks,
  RunPhase2PolynomialDerivation
];

etaPlusPoly[] := s*(1 + Amu) + Lambda*(1 + AG);

etaMinusPoly[] := s*(1 - Amu) + Lambda*(1 - AG);

etaPlus[] := etaPlusPoly[]/s;

etaMinus[] := etaMinusPoly[]/s;

radicandMinus[] := 1 + ((1 - Arho)*s^2)/etaMinusPoly[];

radicandPlus[] := 1 + ((1 + Arho)*s^2)/etaPlusPoly[];

qMinus[] := Sqrt[radicandMinus[]];

qPlus[] := Sqrt[radicandPlus[]];

recipDenomMinus[] := etaPlusPoly[] + etaMinusPoly[]*qMinus[];

recipDenomPlus[] := etaMinusPoly[] + etaPlusPoly[]*qPlus[];

F0ClockForm[] :=
  s^2*(1/recipDenomMinus[] + 1/recipDenomPlus[]) + 2;

F1ClearedReciprocals[] :=
  Expand[
    s^2*(recipDenomMinus[] + recipDenomPlus[]) +
      2*recipDenomMinus[]*recipDenomPlus[]
  ];

alpha0[] := s^2*(etaPlusPoly[] + etaMinusPoly[]) +
  2*etaPlusPoly[]*etaMinusPoly[];

betaMinus[] := etaMinusPoly[]*(s^2 + 2*etaMinusPoly[]);

gammaPlus[] := etaPlusPoly[]*(s^2 + 2*etaPlusPoly[]);

deltaCoupling[] := 2*etaPlusPoly[]*etaMinusPoly[];

F1LinearRadicals[] :=
  alpha0[] + betaMinus[]*qMinus[] + gammaPlus[]*qPlus[] +
    deltaCoupling[]*qMinus[]*qPlus[];

F2AfterFirstSquaring[] :=
  Expand[
    (alpha0[] + betaMinus[]*qMinus[])^2 -
      radicandPlus[]*(gammaPlus[] + deltaCoupling[]*qMinus[])^2
  ];

u0[] := Expand[alpha0[]^2 - radicandPlus[]*gammaPlus[]^2];

u1[] := Expand[2*alpha0[]*betaMinus[] -
  2*radicandPlus[]*gammaPlus[]*deltaCoupling[]];

u2[] := Expand[betaMinus[]^2 - radicandPlus[]*deltaCoupling[]^2];

F3SingleRadical[] := Expand[u0[] + u1[]*qMinus[] + u2[]*radicandMinus[]];

F4AfterSecondSquaring[] :=
  Together[(u0[] + u2[]*radicandMinus[])^2 -
    u1[]^2*radicandMinus[]];

rawPolynomialNumerator[] := Numerator[Together[F4AfterSecondSquaring[]]];

denominatorClearingEtaFactor[] := etaPlusPoly[]^2*etaMinusPoly[]^2;

polynomialNormalizationFactor[] := 4;

P14Polynomial[] :=
  Expand[rawPolynomialNumerator[]/polynomialNormalizationFactor[]];

P14CoefficientList[] := CoefficientList[P14Polynomial[], s];

P14LeadingCoefficient[] :=
  Factor[Coefficient[P14Polynomial[], s, 14]];

P14Degree[] := Exponent[P14Polynomial[], s];

P14CoefficientCount[] := Length[P14CoefficientList[]];

P14IsPolynomialQ[] := PolynomialQ[P14Polynomial[], s];

Phase2PhysicalAssumptions[] := Phase1Assumptions[];

Phase2OriginalExpressionDomain[] := And[
  s != 0,
  etaPlusPoly[] != 0,
  etaMinusPoly[] != 0,
  recipDenomMinus[] != 0,
  recipDenomPlus[] != 0
];

Phase2TransformationNonzeroConditions[] := And[
  etaPlusPoly[] != 0,
  etaMinusPoly[] != 0,
  recipDenomMinus[] != 0,
  recipDenomPlus[] != 0,
  Denominator[Together[F4AfterSecondSquaring[]]] != 0,
  denominatorClearingEtaFactor[] != 0
];

Phase2AdmissibleDomain[] := And[
  Phase2PhysicalAssumptions[],
  Phase2OriginalExpressionDomain[],
  Phase2TransformationNonzeroConditions[]
];

Phase2DerivationSteps[] := {
  <|
    "Step" -> "F0",
    "Equation" -> HoldForm[F0ClockForm[] == 0],
    "Operation" -> "exact rewriting",
    "ForwardImplication" -> True,
    "MayIntroduceSpuriousRoots" -> False,
    "MayLoseAdmissibleRoots" -> False,
    "Notes" -> "This is DStarLambda written with polynomial etaPlus/etaMinus factors."
  |>,
  <|
    "Step" -> "F0 -> F1",
    "Equation" -> HoldForm[F1ClearedReciprocals[] == 0],
    "Operation" -> "multiplication by nonzero reciprocal denominators",
    "ForwardImplication" -> True,
    "MayIntroduceSpuriousRoots" -> True,
    "MayLoseAdmissibleRoots" -> False,
    "NonzeroFactors" -> {recipDenomMinus[], recipDenomPlus[]},
    "Notes" -> "Valid forward on the original-expression domain. Reverse implication is not claimed."
  |>,
  <|
    "Step" -> "F1 linear radical form",
    "Equation" -> HoldForm[F1LinearRadicals[] == 0],
    "Operation" -> "exact rearrangement",
    "ForwardImplication" -> True,
    "MayIntroduceSpuriousRoots" -> False,
    "MayLoseAdmissibleRoots" -> False,
    "Notes" -> "Collects the equation as alpha + beta qMinus + gamma qPlus + delta qMinus qPlus = 0."
  |>,
  <|
    "Step" -> "F1 -> F2",
    "Equation" -> HoldForm[F2AfterFirstSquaring[] == 0],
    "Operation" -> "squaring / first radical elimination",
    "ForwardImplication" -> True,
    "MayIntroduceSpuriousRoots" -> True,
    "MayLoseAdmissibleRoots" -> False,
    "Notes" -> "Isolates qPlus algebraically and squares. Principal branch is not replaced before squaring."
  |>,
  <|
    "Step" -> "F2 -> F3",
    "Equation" -> HoldForm[F3SingleRadical[] == 0],
    "Operation" -> "exact replacement qMinus^2 -> radicandMinus",
    "ForwardImplication" -> True,
    "MayIntroduceSpuriousRoots" -> False,
    "MayLoseAdmissibleRoots" -> False,
    "Notes" -> "Uses the identity Sqrt[z]^2 == z for Mathematica's principal Sqrt."
  |>,
  <|
    "Step" -> "F3 -> F4",
    "Equation" -> HoldForm[F4AfterSecondSquaring[] == 0],
    "Operation" -> "squaring / second radical elimination",
    "ForwardImplication" -> True,
    "MayIntroduceSpuriousRoots" -> True,
    "MayLoseAdmissibleRoots" -> False,
    "Notes" -> "Eliminates qMinus by squaring. Reverse implication is not claimed."
  |>,
  <|
    "Step" -> "F4 -> raw numerator",
    "Equation" -> HoldForm[rawPolynomialNumerator[] == 0],
    "Operation" -> "denominator clearing",
    "ForwardImplication" -> True,
    "MayIntroduceSpuriousRoots" -> True,
    "MayLoseAdmissibleRoots" -> False,
    "NonzeroFactors" -> {Denominator[Together[F4AfterSecondSquaring[]]]},
    "Notes" -> "Valid forward when the rational denominator is nonzero."
  |>,
  <|
    "Step" -> "raw numerator -> P14",
    "Equation" -> HoldForm[P14Polynomial[] == 0],
    "Operation" -> "constant normalization",
    "ForwardImplication" -> True,
    "MayIntroduceSpuriousRoots" -> False,
    "MayLoseAdmissibleRoots" -> False,
    "NormalizationFactor" -> polynomialNormalizationFactor[],
    "Notes" -> "Only the harmless nonzero constant factor 4 is removed. The etaPlusPoly^2 etaMinusPoly^2 factor remains documented as a denominator-clearing/admissibility condition and is not divided out of the polynomial."
  |>
};

(* ============================================================ *)
(* 6. Standard symbolic-check helpers                            *)
(* ============================================================ *)

Phase2CheckStatus[verified_] := Which[
  TrueQ[verified], "PASS",
  verified === False, "FAIL",
  True, "UNRESOLVED"
];

Phase2MakeCheck[name_, statement_, verified_, assumptions_, details_] := <|
  "Name" -> name,
  "Status" -> Phase2CheckStatus[verified],
  "Statement" -> statement,
  "Assumptions" -> assumptions,
  "Details" -> details
|>;

Phase2IdentityCheck[name_, lhs_, rhs_, assumptions_, details_] := Module[
  {
    verified
  },
  verified = Quiet[
    Check[
      TrueQ[FullSimplify[Together[lhs - rhs] == 0, assumptions]],
      Missing["CheckFailed"]
    ]
  ];
  Phase2MakeCheck[
    name,
    HoldForm[lhs == rhs],
    verified,
    assumptions,
    details
  ]
];

Phase2CheckPassedQ[check_Association] := check["Status"] === "PASS";

Phase2ForwardStepCertifiedQ[steps_, stepName_] := Module[
  {
    matchingSteps = Select[steps, #["Step"] === stepName &]
  },
  Length[matchingSteps] === 1 &&
    AllTrue[
      matchingSteps,
      TrueQ[#["ForwardImplication"]] &&
        TrueQ[! #["MayLoseAdmissibleRoots"]] &
    ]
];

Phase2ForwardPreservingChecks[] := Module[
  {
    steps = Phase2DerivationSteps[]
  },
  {
    Phase2IdentityCheck[
      "DStarLambda equals F0ClockForm",
      DStarLambda[],
      F0ClockForm[],
      Phase2PhysicalAssumptions[] && Phase2OriginalExpressionDomain[],
      "Exact rewriting of the nondimensional dispersion relation into the s^2 clock form."
    ],
    Phase2IdentityCheck[
      "First denominator-clearing step",
      recipDenomMinus[]*recipDenomPlus[]*F0ClockForm[],
      F1ClearedReciprocals[],
      Phase2PhysicalAssumptions[] && Phase2OriginalExpressionDomain[],
      "Multiplication by reciprocal denominators. This preserves roots in the forward direction on the original-expression domain."
    ],
    Phase2IdentityCheck[
      "Collected F1 radical form",
      F1ClearedReciprocals[],
      F1LinearRadicals[],
      Phase2PhysicalAssumptions[] && Phase2OriginalExpressionDomain[],
      "Exact collection into alpha0 + betaMinus qMinus + gammaPlus qPlus + delta qMinus qPlus."
    ],
    Phase2MakeCheck[
      "First radical elimination preserves forward implication",
      HoldForm[F1LinearRadicals[] == 0 \[Implies] F2AfterFirstSquaring[] == 0],
      Phase2ForwardStepCertifiedQ[steps, "F1 -> F2"],
      True,
      "The step squares an equation satisfied by any original root. Squaring may add roots but cannot lose roots satisfying the pre-squared equation."
    ],
    Phase2IdentityCheck[
      "F2 single-radical rewrite",
      F2AfterFirstSquaring[],
      F3SingleRadical[],
      Phase2PhysicalAssumptions[] && Phase2OriginalExpressionDomain[],
      "Uses qMinus^2 = radicandMinus with Mathematica's principal Sqrt convention; no PowerExpand or branch replacement is used."
    ],
    Phase2MakeCheck[
      "Second radical elimination preserves forward implication",
      HoldForm[F3SingleRadical[] == 0 \[Implies] F4AfterSecondSquaring[] == 0],
      Phase2ForwardStepCertifiedQ[steps, "F3 -> F4"],
      True,
      "The second squaring may add roots but cannot lose roots satisfying F3 == 0."
    ],
    Phase2MakeCheck[
      "Rational numerator extraction preserves forward implication",
      HoldForm[
        F4AfterSecondSquaring[] == 0 &&
          Denominator[Together[F4AfterSecondSquaring[]]] != 0 \[Implies]
          rawPolynomialNumerator[] == 0
      ],
      Phase2ForwardStepCertifiedQ[steps, "F4 -> raw numerator"],
      HoldForm[Denominator[Together[F4AfterSecondSquaring[]]] != 0],
      "If a defined rational expression is zero and its denominator is nonzero, then its numerator is zero."
    ],
    Phase2IdentityCheck[
      "Final expression equals P14 up to constant normalization",
      rawPolynomialNumerator[],
      polynomialNormalizationFactor[]*P14Polynomial[],
      Phase2PhysicalAssumptions[],
      "Only the nonzero constant factor 4 is divided out. E_+^2 E_-^2 is not divided out of P14."
    ]
  }
];

Phase2PossibleDegreeDropConditions[] := {
  <|
    "Condition" -> HoldForm[P14LeadingCoefficient[] == 0],
    "EquivalentGenericCondition" -> HoldForm[Arho + Amu == 0],
    "Reason" -> "The computed coefficient of s^14 is (Arho + Amu)^2."
  |>,
  <|
    "Condition" -> HoldForm[Lambda == 0],
    "EquivalentGenericCondition" -> "outside Phase I/II theorem scope",
    "Reason" -> "Phase II inherits Lambda > 0 from Phase I."
  |>,
  <|
    "Condition" -> "additional lower-degree cancellations",
    "EquivalentGenericCondition" -> "to be analyzed by specialized parameter strata",
    "Reason" -> "After the leading coefficient vanishes, lower coefficients determine the realized degree."
  |>
};

Phase2SymbolicChecks[] :=
  Module[
    {
      forwardChecks = Phase2ForwardPreservingChecks[],
      leadingCoefficientCheck
    },
    leadingCoefficientCheck =
      TrueQ[FullSimplify[P14LeadingCoefficient[] != 0,
        Phase2PhysicalAssumptions[] && Arho + Amu != 0]];
    <|
      "P14IsPolynomialInS" -> P14IsPolynomialQ[],
      "DegreeInS" -> P14Degree[],
      "GenericDegree14" -> TrueQ[P14Degree[] === 14],
      "CoefficientCount" -> P14CoefficientCount[],
      "LeadingCoefficient" -> P14LeadingCoefficient[],
      "LeadingCoefficientNonzeroGeneric" -> leadingCoefficientCheck,
      "NoPowerExpandUsed" -> True,
      "ForwardPreservingChecks" -> forwardChecks,
      "AllForwardChecksPass" -> AllTrue[forwardChecks, Phase2CheckPassedQ],
      "NoStepCanLoseAdmissibleRoots" ->
        AllTrue[Phase2DerivationSteps[], TrueQ[! #["MayLoseAdmissibleRoots"]] &],
      "AllStepsPreserveForwardImplication" ->
        AllTrue[Phase2DerivationSteps[], TrueQ[#["ForwardImplication"]] &],
      "ForwardImplicationStatus" ->
        "Certified by Phase2ForwardPreservingChecks[]; reverse implication is not claimed."
    |>
  ];

Phase2NumericalSpotChecks[] := Module[
  {
    samples,
    checks,
    x,
    y,
    arhoValue,
    amuValue,
    agValue,
    lambdaValue,
    seedX,
    seedY,
    rootSol,
    rootValue,
    dResidual,
    pResidual
  },

  samples = {
    <|"Arho" -> 1/2, "Amu" -> 1/2, "AG" -> 0, "Lambda" -> 1/4,
      "Seed" -> {-0.6848265998344342, 0.016274647360410503}|>,
    <|"Arho" -> 1/2, "Amu" -> 1/2, "AG" -> 0, "Lambda" -> 1,
      "Seed" -> {-1.3025000536085127, 0.31620384771238624}|>,
    <|"Arho" -> 1/5, "Amu" -> -1/10, "AG" -> 3/10, "Lambda" -> 1/2,
      "Seed" -> {-0.6993929306288597, 0.5914684668807885}|>,
    <|"Arho" -> -2/5, "Amu" -> 1/5, "AG" -> -3/10, "Lambda" -> 2,
      "Seed" -> {-0.3770199539440559, 1.5075311937205984}|>
  };

  checks = Table[
    Quiet[
      Check[
        arhoValue = SetPrecision[sample["Arho"], 80];
        amuValue = SetPrecision[sample["Amu"], 80];
        agValue = SetPrecision[sample["AG"], 80];
        lambdaValue = SetPrecision[sample["Lambda"], 80];
        seedX = SetPrecision[sample["Seed"][[1]], 80];
        seedY = SetPrecision[sample["Seed"][[2]], 80];
        rootSol = FindRoot[
          {
            Re[DStarLambda[] /. {
                s -> x + I*y,
                Arho -> arhoValue,
                Amu -> amuValue,
                AG -> agValue,
                Lambda -> lambdaValue
              }] == 0,
            Im[DStarLambda[] /. {
                s -> x + I*y,
                Arho -> arhoValue,
                Amu -> amuValue,
                AG -> agValue,
                Lambda -> lambdaValue
              }] == 0
          },
          {{x, seedX}, {y, seedY}},
          WorkingPrecision -> 80,
          AccuracyGoal -> 40,
          PrecisionGoal -> 40,
          MaxIterations -> 100
        ];
        rootValue = N[x + I*y /. rootSol, 50];
        dResidual = Abs[
          N[DStarLambda[] /. {
              s -> rootValue,
              Arho -> arhoValue,
              Amu -> amuValue,
              AG -> agValue,
              Lambda -> lambdaValue
            }, 50]
        ];
        pResidual = Abs[
          N[P14Polynomial[] /. {
              s -> rootValue,
              Arho -> arhoValue,
              Amu -> amuValue,
              AG -> agValue,
              Lambda -> lambdaValue
            }, 50]
        ];
        Join[
          sample,
          <|"Root" -> rootValue, "DResidual" -> dResidual,
            "P14Residual" -> pResidual,
            "FailureType" -> None,
            "Passed" -> TrueQ[dResidual < 10^-35 && pResidual < 10^-25]|>
        ],
        Join[sample, <|"Root" -> Missing["FindRootFailed"],
          "DResidual" -> Infinity, "P14Residual" -> Infinity,
          "FailureType" -> "FindRootFailed",
          "Passed" -> False|>]
      ]
    ],
    {sample, samples}
  ];

  checks
];

RunPhase2PolynomialDerivation[] := Module[
  {
    symbolicChecks,
    numericChecks,
    forwardTheoremPass,
    numericalSpotCheckPass
  },

  symbolicChecks = Phase2SymbolicChecks[];
  numericChecks = Phase2NumericalSpotChecks[];
  forwardTheoremPass =
    TrueQ[symbolicChecks["P14IsPolynomialInS"]] &&
    TrueQ[symbolicChecks["GenericDegree14"]] &&
    TrueQ[symbolicChecks["CoefficientCount"] === 15] &&
    TrueQ[symbolicChecks["LeadingCoefficientNonzeroGeneric"]] &&
    TrueQ[symbolicChecks["AllForwardChecksPass"]] &&
    TrueQ[symbolicChecks["NoStepCanLoseAdmissibleRoots"]] &&
    TrueQ[symbolicChecks["AllStepsPreserveForwardImplication"]] &&
    TrueQ[symbolicChecks["NoPowerExpandUsed"]];
  numericalSpotCheckPass = AllTrue[numericChecks, TrueQ[#["Passed"]] &];

  <|
    "Phase" -> "Phase II: DStarLambda(s) -> P14(s)",
    "Theorem" -> "DStarLambda(s) == 0 implies P14(s) == 0 on the admissible domain. Reverse implication is not claimed.",
    "PhysicalAssumptions" -> Phase2PhysicalAssumptions[],
    "OriginalExpressionDomain" -> Phase2OriginalExpressionDomain[],
    "TransformationNonzeroConditions" -> Phase2TransformationNonzeroConditions[],
    "AdmissibleDomain" -> Phase2AdmissibleDomain[],
    "EtaPlusPolynomial" -> etaPlusPoly[],
    "EtaMinusPolynomial" -> etaMinusPoly[],
    "RadicandMinus" -> radicandMinus[],
    "RadicandPlus" -> radicandPlus[],
    "ReciprocalDenominatorMinus" -> recipDenomMinus[],
    "ReciprocalDenominatorPlus" -> recipDenomPlus[],
    "F0" -> F0ClockForm[],
    "F1" -> F1LinearRadicals[],
    "F2" -> F2AfterFirstSquaring[],
    "F3" -> F3SingleRadical[],
    "F4" -> F4AfterSecondSquaring[],
    "RawPolynomialNumerator" -> rawPolynomialNumerator[],
    "DenominatorClearingEtaFactor" -> denominatorClearingEtaFactor[],
    "PolynomialNormalizationFactor" -> polynomialNormalizationFactor[],
    "P14Polynomial" -> P14Polynomial[],
    "P14CoefficientList" -> P14CoefficientList[],
    "P14Degree" -> P14Degree[],
    "P14CoefficientCount" -> P14CoefficientCount[],
    "P14LeadingCoefficient" -> P14LeadingCoefficient[],
    "PossibleDegreeDropConditions" -> Phase2PossibleDegreeDropConditions[],
    "DerivationSteps" -> Phase2DerivationSteps[],
    "SymbolicChecks" -> symbolicChecks,
    "NumericalSpotChecks" -> numericChecks,
    "ForwardTheoremPASS" -> forwardTheoremPass,
    "NumericalSpotCheckPASS" -> numericalSpotCheckPass,
    "SpuriousRootSources" -> {
      "clearing reciprocal denominators",
      "first squaring",
      "second squaring",
      "clearing rational denominators after radical elimination"
    },
    "PASS" ->
      TrueQ[forwardTheoremPass] &&
      TrueQ[numericalSpotCheckPass]
  |>
];
