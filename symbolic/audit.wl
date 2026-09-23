(* ::Package:: *)

(* ============================================================ *)
(* Phase I audit entry point                                    *)
(* ============================================================ *)

(* Purpose:
   Run the frozen Phase I dimensional-to-nondimensional verification audit.

   Command-line use:

     wolframscript -script symbolic/audit.wl

   Fresh-notebook use:

     Get[FileNameJoin[{projectRoot, "symbolic", "audit.wl"}]]

   In a notebook, this file returns an Association and does not call Exit[].
   In wolframscript, it prints a PASS/FAIL summary and exits with the
   corresponding status code.
*)

ClearAll[
  SymbolicAuditDirectory,
  RunSymbolicAudit,
  SymbolicAuditResult,
  SymbolicAuditCommandLineQ
];

SymbolicAuditDirectory[] :=
  If[StringQ[$InputFileName] && $InputFileName =!= "",
    DirectoryName[$InputFileName],
    If[ValueQ[projectRoot],
      FileNameJoin[{projectRoot, "symbolic"}],
      Directory[]
    ]
  ];

RunSymbolicAudit[] := Module[
  {
    scriptDirectory,
    logDirectory,
    exportDirectory,
    timestamp,
    result,
    allChecksPassed,
    logPath,
    resultPath,
    dStarTeXPath,
    scaledTeXPath
  },

  scriptDirectory = SymbolicAuditDirectory[];

  Get[FileNameJoin[{scriptDirectory, "verification", "Phase1Verification.wl"}]];

  logDirectory = FileNameJoin[{scriptDirectory, "logs"}];
  exportDirectory = FileNameJoin[{scriptDirectory, "exports"}];
  If[! DirectoryQ[logDirectory],
    CreateDirectory[logDirectory, CreateIntermediateDirectories -> True]
  ];
  If[! DirectoryQ[exportDirectory],
    CreateDirectory[exportDirectory, CreateIntermediateDirectories -> True]
  ];

  timestamp =
    DateString[{"Year", "Month", "Day", "_", "Hour", "Minute", "Second"}];

  result = RunPhase1Verification[];
  allChecksPassed =
    TrueQ[result["EquivalenceVerified"]] &&
      AllTrue[result["VerificationTable"], TrueQ[#["Passed"]] &];

  logPath = FileNameJoin[{logDirectory, "phase1_audit_" <> timestamp <> ".txt"}];
  resultPath = FileNameJoin[{logDirectory, "phase1_result_" <> timestamp <> ".wl"}];
  dStarTeXPath =
    FileNameJoin[{exportDirectory, "DStarLambda_phase1_" <> timestamp <> ".tex"}];
  scaledTeXPath =
    FileNameJoin[
      {exportDirectory, "ScaledDimensionalD_phase1_" <> timestamp <> ".tex"}
    ];

  Export[
    logPath,
    StringRiffle[
      {
        "Phase I symbolic dispersion-relation audit",
        "Timestamp: " <> DateString[],
        "Verification target: DStarLambda == rhoT/(2 k^2) * D(gamma)",
        "Equivalence verified: " <> ToString[result["EquivalenceVerified"]],
        "All component checks passed: " <> ToString[allChecksPassed],
        "Scaling prefactor: " <> ToString[result["ScalingPrefactor"], InputForm],
        "Eta scaling differences: " <> ToString[result["EtaScalingDifferences"], InputForm],
        "Radicand mapping differences: " <> ToString[result["RadicandMappingDifferences"], InputForm],
        "Radical mapping differences: " <> ToString[result["RadicalMappingDifferences"], InputForm],
        "Denominator scaling differences: " <> ToString[result["DenominatorScalingDifferences"], InputForm],
        "Constant term difference: " <> ToString[result["ConstantTermDifference"], InputForm],
        "Reciprocal term difference: " <> ToString[result["ReciprocalTermDifference"], InputForm],
        "Difference: " <> ToString[result["Difference"], InputForm],
        "Assumptions: " <> ToString[result["Assumptions"], InputForm],
        "Not yet performed: " <> ToString[result["NotYetPerformed"], InputForm]
      },
      "\n"
    ],
    "Text"
  ];

  Put[result, resultPath];

  Export[
    dStarTeXPath,
    ToString[TeXForm[result["DStarLambda"]]],
    "Text"
  ];

  Export[
    scaledTeXPath,
    ToString[TeXForm[result["ScaledDimensionalD"]]],
    "Text"
  ];

  Join[
    result,
    <|
      "LogPath" -> logPath,
      "ResultPath" -> resultPath,
      "DStarTeXPath" -> dStarTeXPath,
      "ScaledDimensionalDTeXPath" -> scaledTeXPath,
      "PassFailSummary" ->
        If[allChecksPassed,
          "PASS: Phase I symbolic nondimensionalization verified.",
          "FAIL: at least one Phase I symbolic check did not verify."
        ]
    |>
  ]
];

SymbolicAuditResult = RunSymbolicAudit[];

SymbolicAuditCommandLineQ := TrueQ[$FrontEnd === Null];

If[SymbolicAuditCommandLineQ,
  Print[SymbolicAuditResult["PassFailSummary"]];
  Print["Log: ", SymbolicAuditResult["LogPath"]];
  If[TrueQ[SymbolicAuditResult["EquivalenceVerified"]],
    Exit[0],
    Print["Difference: ", SymbolicAuditResult["Difference"]];
    Exit[1]
  ],
  SymbolicAuditResult
]
