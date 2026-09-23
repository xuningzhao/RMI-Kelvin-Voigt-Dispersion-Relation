(* ::Package:: *)

(* ============================================================ *)
(* Phase III audit entry point                                  *)
(* ============================================================ *)

(* Purpose:
   Run the Phase III polynomial-completeness/no-root-loss audit.

     OmegaDefinition && DStarLambda(s) == 0  ==>  P14Polynomial(s) == 0

   Dependencies:
     - verification/Phase1Verification.wl
     - verification/Phase2PolynomialDerivation.wl
     - verification/Phase3PolynomialCompleteness.wl

   Expected behavior:
     - PASS means the Phase II forward theorem and supplemental numerical
       validation both passed.
     - FAIL includes unresolved symbolic checks.
     - The reverse implication is not claimed.

   Candidate filtering and physical admissibility are later-stage tasks.
*)

ClearAll[
  SymbolicPhase3AuditDirectory,
  SymbolicPhase3PrintSummary,
  RunSymbolicPhase3Audit,
  SymbolicPhase3AuditResult,
  SymbolicPhase3AuditCommandLineQ
];

If[! ValueQ[$AuditVerbosity], $AuditVerbosity = "Summary"];

SymbolicPhase3AuditDirectory[] :=
  If[StringQ[$InputFileName] && $InputFileName =!= "",
    DirectoryName[$InputFileName],
    If[ValueQ[projectRoot],
      FileNameJoin[{projectRoot, "symbolic"}],
      Directory[]
    ]
  ];

SymbolicPhase3PrintSummary[result_Association] := Module[
  {
    checks = result["Phase2ForwardPreservingChecks"]
  },
  If[$AuditVerbosity === "Quiet", Return[Null]];
  Scan[
    Print["[" <> #["Status"] <> "] " <> #["Name"]] &,
    checks
  ];
  Print[
    "Overall Phase III status: ",
    If[TrueQ[result["PASS"]], "PASS", "FAIL"]
  ];
];

RunSymbolicPhase3Audit[] := Module[
  {
    scriptDirectory,
    logDirectory,
    exportDirectory,
    timestamp,
    result,
    reportPath,
    logPath,
    resultPath,
    noRootLossPath,
    forwardChecksPath,
    numericalValidationPath,
    statusDisplay
  },

  scriptDirectory = SymbolicPhase3AuditDirectory[];

  Get[FileNameJoin[{scriptDirectory, "verification", "Phase1Verification.wl"}]];
  Get[FileNameJoin[{scriptDirectory, "verification", "Phase2PolynomialDerivation.wl"}]];
  Get[FileNameJoin[{scriptDirectory, "verification", "Phase3PolynomialCompleteness.wl"}]];

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

  result = RunPhase3PolynomialCompleteness[];
  SymbolicPhase3PrintSummary[result];

  statusDisplay =
    Grid[
      {
        {"Forward implication:",
          Style[
            If[TrueQ[result["ForwardImplicationPASS"]], "PASS", "FAIL"],
            If[TrueQ[result["ForwardImplicationPASS"]], Darker[Green], Red],
            Bold
          ]},
        {"Numerical validation:",
          Style[
            If[TrueQ[result["NumericalValidationPASS"]], "PASS", "FAIL"],
            If[TrueQ[result["NumericalValidationPASS"]], Darker[Green], Red],
            Bold
          ]},
        {"Reverse implication claimed:",
          Style[
            If[TrueQ[result["ReverseImplicationClaimed"]], "YES", "NO"],
            If[TrueQ[result["ReverseImplicationClaimed"]], Red, Darker[Green]],
            Bold
          ]},
        {"Overall Phase III:",
          Style[
            result["Phase3Status"],
            If[result["Phase3Status"] === "PASS", Darker[Green], Red],
            Bold,
            16
          ]}
      },
      Alignment -> Left,
      Spacings -> {2, 0.8}
    ];

  result = Join[result, <|"Phase3StatusDisplay" -> statusDisplay|>];

  reportPath =
    FileNameJoin[{exportDirectory,
      "phase3_polynomial_completeness_report_" <> timestamp <> ".md"}];
  logPath =
    FileNameJoin[{logDirectory, "phase3_audit_" <> timestamp <> ".txt"}];
  resultPath =
    FileNameJoin[{logDirectory, "phase3_result_" <> timestamp <> ".wl"}];
  noRootLossPath =
    FileNameJoin[{exportDirectory,
      "NoRootLossConditions_" <> timestamp <> ".wl"}];
  forwardChecksPath =
    FileNameJoin[{exportDirectory,
      "Phase2ForwardCompletenessChecks_" <> timestamp <> ".wl"}];
  numericalValidationPath =
    FileNameJoin[{exportDirectory,
      "phase3_forward_numerical_validation_" <> timestamp <> ".wl"}];

  Export[reportPath, Phase3CompletenessReportString[result], "Text"];
  Put[result, resultPath];
  Put[result["NoRootLossConditions"], noRootLossPath];
  Put[result["Phase2ForwardPreservingChecks"], forwardChecksPath];
  Put[result["NumericalForwardValidation"], numericalValidationPath];

  Export[
    logPath,
    StringRiffle[
      {
        "Phase III polynomial-completeness audit",
        "Timestamp: " <> DateString[],
        "PASS: " <> ToString[result["PASS"]],
        "Phase3Status: " <> ToString[result["Phase3Status"]],
        "ForwardImplicationPASS: " <> ToString[result["ForwardImplicationPASS"]],
        "PolynomialCompletenessPASS: " <> ToString[result["PolynomialCompletenessPASS"]],
        "NumericalValidationPASS: " <> ToString[result["NumericalValidationPASS"]],
        "Reverse implication claimed: " <> ToString[result["ReverseImplicationClaimed"]],
        "Extraneous polynomial roots expected: " <> ToString[result["ExtraneousPolynomialRootsExpected"]],
        "Report: " <> reportPath,
        "NoRootLossConditions export: " <> noRootLossPath,
        "Phase II forward-completeness checks export: " <> forwardChecksPath,
        "Numerical validation export: " <> numericalValidationPath,
        "Result: " <> resultPath,
        "Summary: " <> ToString[result["PassFailSummary"]]
      },
      "\n"
    ],
    "Text"
  ];

  Join[
    result,
    <|
      "ReportPath" -> reportPath,
      "NoRootLossConditionsPath" -> noRootLossPath,
      "Phase2ForwardCompletenessChecksPath" -> forwardChecksPath,
      "NumericalValidationPath" -> numericalValidationPath,
      "LogPath" -> logPath,
      "ResultPath" -> resultPath
    |>
  ]
];

SymbolicPhase3AuditResult = RunSymbolicPhase3Audit[];

SymbolicPhase3AuditCommandLineQ := TrueQ[$FrontEnd === Null];

If[SymbolicPhase3AuditCommandLineQ,
  Print[SymbolicPhase3AuditResult["PassFailSummary"]];
  Print["Log: ", SymbolicPhase3AuditResult["LogPath"]];
  If[TrueQ[SymbolicPhase3AuditResult["PASS"]],
    Exit[0],
    Exit[1]
  ],
  SymbolicPhase3AuditResult
]
