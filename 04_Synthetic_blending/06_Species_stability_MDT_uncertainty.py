
# ============================================================
# ANALYSIS 2
# MDT SAMPLING UNCERTAINTY AND ROBUSTNESS
# CRC MICROBIOME — FIXED SVM OBSERVER
# ============================================================

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from sklearn.pipeline import Pipeline
from sklearn.impute import SimpleImputer
from sklearn.svm import SVC
from pathlib import Path

# ============================================================
# 1. CONFIGURATION
# ============================================================

OUT = Path("MDT_Analysis2")
OUT.mkdir(exist_ok=True)

N_RUNS = 100
N_PAIRS = 50
SEED = 2026

ALPHAS = np.round(np.linspace(0, 1, 21), 8)
DOMINANCE = 0.80

TRANSITIONS = [
    ("Healthy", "Adenoma"),
    ("Healthy", "Cancer"),
    ("Adenoma", "Cancer")
]

STRATEGIES = [
    "fixed_pairs",
    "resample_each_alpha"
]

PRIMARY_MDT = {
    "Healthy -> Adenoma": 0.65,
    "Healthy -> Cancer": 0.60,
    "Adenoma -> Cancer": 0.55
}

# ============================================================
# 2. LOAD AND VALIDATE DATA
# ============================================================

X = pd.read_csv(
    "/content/TPM_clr_batch_corrected_v2.csv",
    index_col=0
)

meta = pd.read_csv(
    "/content/CRC_metadata_aligned_bact.csv"
)

X.index = X.index.astype(str)
meta["Sample_ID"] = meta["Sample_ID"].astype(str)

assert X.index.is_unique
assert meta["Sample_ID"].is_unique

common = X.index.intersection(meta["Sample_ID"])

meta = meta.loc[
    meta["Sample_ID"].isin(common)
].copy()

X = X.loc[meta["Sample_ID"]]

assert np.array_equal(
    X.index.to_numpy(),
    meta["Sample_ID"].to_numpy()
)

X = X.apply(pd.to_numeric, errors="raise")

assert np.isfinite(X.to_numpy()).all()

y = meta["Group"].to_numpy()
features = X.columns.tolist()

print("Samples:", len(X))
print("Features:", X.shape[1])
print(pd.Series(y).value_counts())

assert len(X) == 995
assert X.shape[1] == 1023

# ============================================================
# 3. TRAIN FIXED SVM OBSERVER
# Same configuration as original blending script
# ============================================================

observer = Pipeline([
    ("imputer", SimpleImputer(strategy="median")),
    ("svm", SVC(
        kernel="rbf",
        C=10,
        gamma="scale",
        probability=True,
        class_weight="balanced",
        random_state=42
    ))
])

observer.fit(X, y)

classes = observer.named_steps["svm"].classes_

by_group = {
    group: X.loc[y == group].to_numpy()
    for group in classes
}

# Preserve feature names during prediction
def predict_labels(matrix):
    df = pd.DataFrame(matrix, columns=features)
    return observer.predict(df)

# ============================================================
# 4. ONE SIMULATION RUN
# ============================================================

def simulate_transition(recipient, donor, strategy, seed):

    rng = np.random.default_rng(seed)

    Xr = by_group[recipient]
    Xd = by_group[donor]

    if strategy == "fixed_pairs":
        ri = rng.integers(len(Xr), size=N_PAIRS)
        di = rng.integers(len(Xd), size=N_PAIRS)

        R_fixed = Xr[ri]
        D_fixed = Xd[di]

    records = []
    sampled_distances = []

    for alpha in ALPHAS:

        if strategy == "fixed_pairs":
            R = R_fixed
            D = D_fixed
        else:
            ri = rng.integers(len(Xr), size=N_PAIRS)
            di = rng.integers(len(Xd), size=N_PAIRS)
            R = Xr[ri]
            D = Xd[di]

        blends = (1 - alpha) * R + alpha * D

        labels = predict_labels(blends)

        fraction = np.mean(labels == donor)

        distances = np.linalg.norm(D - R, axis=1)

        sampled_distances.append(distances)

        records.append({
            "Alpha": float(alpha),
            "Donor_fraction": float(fraction),
            "N_donor_labels": int(np.sum(labels == donor)),
            "Endpoint_distance_mean": float(distances.mean()),
            "Endpoint_distance_median": float(np.median(distances)),
            "Recipient_correct_fraction":
                float(np.mean(labels == recipient))
                if np.isclose(alpha, 0) else np.nan,
            "Donor_correct_fraction":
                float(np.mean(labels == donor))
                if np.isclose(alpha, 1) else np.nan
        })

    curve = pd.DataFrame(records)

    hits = curve.loc[
        curve["Donor_fraction"] >= DOMINANCE
    ]

    if len(hits):
        threshold_row = hits.iloc[0]
        mdt = float(threshold_row["Alpha"])

        # Distance of paths sampled at the actual threshold
        endpoint_mean = float(
            threshold_row["Endpoint_distance_mean"]
        )

        displacement = mdt * endpoint_mean
    else:
        mdt = np.nan
        endpoint_mean = np.nan
        displacement = np.nan

    # Overall sampled path-distance diagnostic
    all_distances = np.concatenate(sampled_distances)

    summary = {
        "MDT": mdt,
        "Endpoint_distance_at_MDT": endpoint_mean,
        "Displacement_at_MDT": displacement,
        "All_sampled_endpoint_distance_mean":
            float(all_distances.mean()),
        "Donor_fraction_alpha0":
            float(curve.iloc[0]["Donor_fraction"]),
        "Donor_fraction_alpha1":
            float(curve.iloc[-1]["Donor_fraction"]),
        "Recipient_correct_alpha0":
            float(curve.iloc[0]["Recipient_correct_fraction"]),
        "Donor_correct_alpha1":
            float(curve.iloc[-1]["Donor_correct_fraction"])
    }

    return summary, curve

# ============================================================
# 5. RUN REPEATED SIMULATIONS
# ============================================================

run_rows = []
curve_rows = []

# Same seed index across transitions enables paired
# comparison of simulation runs. Samples remain independent.
for strategy_idx, strategy in enumerate(STRATEGIES):

    print("\nStrategy:", strategy)

    for run in range(N_RUNS):

        for t_idx, (recipient, donor) in enumerate(TRANSITIONS):

            seed = (
                SEED
                + strategy_idx * 1000000
                + run * 100
                + t_idx
            )

            result, curve = simulate_transition(
                recipient, donor, strategy, seed
            )

            transition = f"{recipient} -> {donor}"

            run_rows.append({
                "Strategy": strategy,
                "Run": run + 1,
                "Seed": seed,
                "Transition": transition,
                "Original_primary_MDT":
                    PRIMARY_MDT[transition],
                **result
            })

            curve["Strategy"] = strategy
            curve["Run"] = run + 1
            curve["Transition"] = transition

            curve_rows.append(curve)

        if (run + 1) % 10 == 0:
            print(f"Completed {run + 1}/{N_RUNS}")

runs = pd.DataFrame(run_rows)
curves = pd.concat(curve_rows, ignore_index=True)

runs.to_csv(
    OUT / "MDT_100Run_Results.csv",
    index=False
)

curves.to_csv(
    OUT / "MDT_100Run_Dominance_Curves.csv",
    index=False
)

# ============================================================
# 6. SUMMARIZE MDT AND DISPLACEMENT UNCERTAINTY
# ============================================================

def summarize(x):
    x = np.asarray(x, dtype=float)
    x = x[np.isfinite(x)]

    if not len(x):
        return pd.Series({
            "Mean": np.nan,
            "Median": np.nan,
            "SD": np.nan,
            "CI_low": np.nan,
            "CI_high": np.nan,
            "Min": np.nan,
            "Max": np.nan,
            "N_valid": 0
        })

    return pd.Series({
        "Mean": np.mean(x),
        "Median": np.median(x),
        "SD": np.std(x, ddof=1) if len(x) > 1 else np.nan,
        "CI_low": np.quantile(x, 0.025),
        "CI_high": np.quantile(x, 0.975),
        "Min": np.min(x),
        "Max": np.max(x),
        "N_valid": len(x)
    })

mdt_summary = (
    runs.groupby(["Strategy", "Transition"])["MDT"]
    .apply(summarize)
    .unstack()
    .reset_index()
)

distance_summary = (
    runs.groupby(
        ["Strategy", "Transition"]
    )["Displacement_at_MDT"]
    .apply(summarize)
    .unstack()
    .reset_index()
)

mdt_summary.to_csv(
    OUT / "MDT_Uncertainty_Summary.csv",
    index=False
)

distance_summary.to_csv(
    OUT / "MDT_Displacement_Uncertainty.csv",
    index=False
)

print("\n=== MDT UNCERTAINTY ===")
print(mdt_summary.round(4).to_string(index=False))

print("\n=== DISPLACEMENT UNCERTAINTY ===")
print(distance_summary.round(4).to_string(index=False))

# ============================================================
# 7. PAIRED TRANSITION COMPARISONS
# ============================================================

comparison_rows = []
ordering_rows = []

for strategy in STRATEGIES:

    sub = runs.loc[runs["Strategy"] == strategy]

    wide = sub.pivot(
        index="Run",
        columns="Transition",
        values="MDT"
    )

    A = wide["Healthy -> Adenoma"]
    B = wide["Healthy -> Cancer"]
    C = wide["Adenoma -> Cancer"]

    comparisons = [
        ("Healthy -> Adenoma", "Healthy -> Cancer", A, B),
        ("Healthy -> Adenoma", "Adenoma -> Cancer", A, C),
        ("Healthy -> Cancer", "Adenoma -> Cancer", B, C)
    ]

    for name1, name2, x, z in comparisons:

        difference = (x - z).dropna()

        comparison_rows.append({
            "Strategy": strategy,
            "Transition_1": name1,
            "Transition_2": name2,
            "Mean_difference": difference.mean(),
            "Median_difference": difference.median(),
            "CI_low": difference.quantile(0.025),
            "CI_high": difference.quantile(0.975),
            "Fraction_first_lower": np.mean(difference < 0),
            "Fraction_equal": np.mean(difference == 0),
            "Fraction_first_higher": np.mean(difference > 0),
            "N_valid": len(difference)
        })

    valid = wide.dropna()

    if len(valid):
        strict_lowest = (
            (valid["Adenoma -> Cancer"] <
             valid["Healthy -> Adenoma"])
            &
            (valid["Adenoma -> Cancer"] <
             valid["Healthy -> Cancer"])
        )

        tied_lowest = (
            valid["Adenoma -> Cancer"]
            == valid.min(axis=1)
        )

        ordering_rows.append({
            "Strategy": strategy,
            "N_valid_runs": len(valid),
            "Fraction_AC_strictly_lowest":
                strict_lowest.mean(),
            "Fraction_AC_lowest_including_ties":
                tied_lowest.mean(),
            "Fraction_all_three_equal":
                (
                    valid.nunique(axis=1) == 1
                ).mean()
        })

comparisons_df = pd.DataFrame(comparison_rows)
ordering_df = pd.DataFrame(ordering_rows)

comparisons_df.to_csv(
    OUT / "MDT_Pairwise_Differences.csv",
    index=False
)

ordering_df.to_csv(
    OUT / "MDT_Ordering_Stability.csv",
    index=False
)

print("\n=== PAIRWISE MDT DIFFERENCES ===")
print(comparisons_df.round(4).to_string(index=False))

print("\n=== ORDERING STABILITY ===")
print(ordering_df.round(4).to_string(index=False))

# ============================================================
# 8. ENDPOINT CLASSIFICATION DIAGNOSTICS
# ============================================================

endpoint_summary = (
    runs.groupby(["Strategy", "Transition"])
    .agg(
        Recipient_correct_mean=(
            "Recipient_correct_alpha0", "mean"
        ),
        Donor_correct_mean=(
            "Donor_correct_alpha1", "mean"
        ),
        Donor_fraction_alpha0_mean=(
            "Donor_fraction_alpha0", "mean"
        ),
        Donor_fraction_alpha1_mean=(
            "Donor_fraction_alpha1", "mean"
        )
    )
    .reset_index()
)

endpoint_summary.to_csv(
    OUT / "MDT_Endpoint_Classification.csv",
    index=False
)

print("\n=== ENDPOINT DIAGNOSTICS ===")
print(endpoint_summary.round(4).to_string(index=False))

# ============================================================
# 9. FIGURE 1: MDT DISTRIBUTIONS
# ============================================================

colors = {
    "Healthy -> Adenoma": "#984ea3",
    "Healthy -> Cancer": "#4daf4a",
    "Adenoma -> Cancer": "#d4af37"
}

for strategy in STRATEGIES:

    fig, ax = plt.subplots(figsize=(8, 5))

    values = [
        runs.loc[
            (runs.Strategy == strategy) &
            (runs.Transition == transition),
            "MDT"
        ].dropna().values
        for transition in PRIMARY_MDT
    ]

    bp = ax.boxplot(
        values,
        tick_labels=list(PRIMARY_MDT),
        patch_artist=True,
        widths=0.5
    )

    for patch, transition in zip(
        bp["boxes"], PRIMARY_MDT
    ):
        patch.set_facecolor(colors[transition])
        patch.set_alpha(0.7)

    ax.set_ylabel("Minimal dominance threshold (alpha)")
    ax.set_title(f"MDT sampling uncertainty: {strategy}")
    ax.set_ylim(0, 1)

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    plt.xticks(rotation=15, ha="right")
    plt.tight_layout()

    plt.savefig(
        OUT / f"MDT_Distribution_{strategy}.png",
        dpi=600,
        bbox_inches="tight"
    )
    plt.show()

# ============================================================
# 10. FIGURE 2: DISTANCE-ADJUSTED DISPLACEMENT
# ============================================================

for strategy in STRATEGIES:

    fig, ax = plt.subplots(figsize=(8, 5))

    values = [
        runs.loc[
            (runs.Strategy == strategy) &
            (runs.Transition == transition),
            "Displacement_at_MDT"
        ].dropna().values
        for transition in PRIMARY_MDT
    ]

    bp = ax.boxplot(
        values,
        tick_labels=list(PRIMARY_MDT),
        patch_artist=True,
        widths=0.5
    )

    for patch, transition in zip(
        bp["boxes"], PRIMARY_MDT
    ):
        patch.set_facecolor(colors[transition])
        patch.set_alpha(0.7)

    ax.set_ylabel(
        "Displacement at MDT\n"
        "(corrected CLR feature-space distance)"
    )
    ax.set_title(
        f"Distance-adjusted thresholds: {strategy}"
    )

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    plt.xticks(rotation=15, ha="right")
    plt.tight_layout()

    plt.savefig(
        OUT / f"MDT_Displacement_{strategy}.png",
        dpi=600,
        bbox_inches="tight"
    )
    plt.show()

# ============================================================
# 11. DOMINANCE CURVE UNCERTAINTY
# ============================================================

curve_summary = (
    curves.groupby(
        ["Strategy", "Transition", "Alpha"]
    )["Donor_fraction"]
    .agg(
        Mean="mean",
        SD="std",
        CI_low=lambda x: x.quantile(0.025),
        CI_high=lambda x: x.quantile(0.975)
    )
    .reset_index()
)

curve_summary.to_csv(
    OUT / "MDT_Dominance_Curve_Uncertainty.csv",
    index=False
)

for strategy in STRATEGIES:

    fig, ax = plt.subplots(figsize=(8, 5))

    for transition in PRIMARY_MDT:

        sub = curve_summary.loc[
            (curve_summary.Strategy == strategy) &
            (curve_summary.Transition == transition)
        ].sort_values("Alpha")

        ax.plot(
            sub["Alpha"].to_numpy(),
            sub["Mean"].to_numpy(),
            label=transition,
            color=colors[transition],
            linewidth=2.5
        )

        ax.fill_between(
            sub["Alpha"].to_numpy(),
            sub["CI_low"].to_numpy(),
            sub["CI_high"].to_numpy(),
            color=colors[transition],
            alpha=0.15
        )

    ax.axhline(
        DOMINANCE,
        linestyle="--",
        color="black",
        linewidth=1
    )

    ax.set_xlabel("Donor mixing fraction (alpha)")
    ax.set_ylabel("Fraction classified as donor")
    ax.set_title(
        f"Donor dominance uncertainty: {strategy}"
    )

    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1.02)
    ax.legend(frameon=False)

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    plt.tight_layout()

    plt.savefig(
        OUT / f"MDT_Dominance_Uncertainty_{strategy}.png",
        dpi=600,
        bbox_inches="tight"
    )
    plt.show()

# ============================================================
# 12. FINAL CHECKS
# ============================================================

assert len(runs) == (
    len(STRATEGIES) * N_RUNS * len(TRANSITIONS)
)

assert runs["MDT"].dropna().between(0, 1).all()

print("\n==========================================")
print("ANALYSIS 2 COMPLETED")
print("==========================================")

print("Simulation runs:", len(runs))
print("Output folder:", OUT.resolve())

print("\nMain files:")
print("- MDT_Uncertainty_Summary.csv")
print("- MDT_Displacement_Uncertainty.csv")
print("- MDT_Pairwise_Differences.csv")
print("- MDT_Ordering_Stability.csv")
print("- MDT_Endpoint_Classification.csv")
print("- MDT_Dominance_Curve_Uncertainty.csv")


# ============================================================
# ANALYSIS 3
# SPECIES PERTURBATION RANK STABILITY
# CRC MDT-TM FRAMEWORK
# ============================================================

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from pathlib import Path
from sklearn.pipeline import Pipeline
from sklearn.impute import SimpleImputer
from sklearn.svm import SVC
from scipy.stats import spearmanr

# ============================================================
# 1. CONFIGURATION
# ============================================================

OUT = Path("MDT_Analysis3_Rank_Stability")
OUT.mkdir(exist_ok=True)

SEED = 2026
N_RUNS = 30             # Change to 100 for final analysis
N_PAIRS = 50
DELTA = 0.05
EPSILON = 0.05
TOP_K = 10

# Number of species perturbed in each prediction batch
FEATURE_CHUNK = 20

TRANSITIONS = [
    ("Healthy", "Adenoma"),
    ("Healthy", "Cancer"),
    ("Adenoma", "Cancer")
]

PRIMARY_MDT = {
    ("Healthy", "Adenoma"): 0.65,
    ("Healthy", "Cancer"): 0.60,
    ("Adenoma", "Cancer"): 0.55
}

# ============================================================
# 2. LOAD AND ALIGN DATA
# ============================================================

X = pd.read_csv(
    "/content/TPM_clr_batch_corrected_v2.csv",
    index_col=0
)

meta = pd.read_csv(
    "/content/CRC_metadata_aligned_bact.csv"
)

X.index = X.index.astype(str)
meta["Sample_ID"] = meta["Sample_ID"].astype(str)

assert X.index.is_unique
assert meta["Sample_ID"].is_unique

common = X.index.intersection(meta["Sample_ID"])

meta = meta[
    meta["Sample_ID"].isin(common)
].copy()

X = X.loc[meta["Sample_ID"]]

assert np.array_equal(
    X.index.to_numpy(),
    meta["Sample_ID"].to_numpy()
)

X = X.apply(pd.to_numeric, errors="raise")

assert np.isfinite(X.to_numpy()).all()

y = meta["Group"].to_numpy()
features = X.columns.to_numpy()

print("Samples:", len(X))
print("Features:", len(features))
print(pd.Series(y).value_counts())

assert len(X) == 995
assert X.shape[1] == 1023

# ============================================================
# 3. FIT ORIGINAL SVM OBSERVER
# ============================================================

observer = Pipeline([
    ("imputer", SimpleImputer(strategy="median")),
    ("svm", SVC(
        kernel="rbf",
        C=10,
        gamma="scale",
        probability=True,
        class_weight="balanced",
        random_state=42
    ))
])

observer.fit(X, y)

classes = observer.named_steps["svm"].classes_

by_group = {
    group: X.loc[y == group].to_numpy()
    for group in classes
}

# Use DataFrames to preserve feature names
def donor_probabilities(matrix, donor_idx):
    df = pd.DataFrame(matrix, columns=features)
    return observer.predict_proba(df)[:, donor_idx]

# ============================================================
# 4. DIRECTIONAL PERTURBATION FUNCTION
# ============================================================

def calculate_influence(recipient, donor, seed):

    rng = np.random.default_rng(seed)

    R_pool = by_group[recipient]
    D_pool = by_group[donor]

    ri = rng.integers(len(R_pool), size=N_PAIRS)
    di = rng.integers(len(D_pool), size=N_PAIRS)

    R = R_pool[ri]
    D = D_pool[di]

    alpha = max(
        0.0,
        PRIMARY_MDT[(recipient, donor)] - EPSILON
    )

    blends = (1 - alpha) * R + alpha * D

    donor_idx = np.where(classes == donor)[0][0]

    n_features = blends.shape[1]
    effects = np.zeros(n_features)

    # Process species in chunks to avoid excessive memory use
    for start in range(0, n_features, FEATURE_CHUNK):

        stop = min(start + FEATURE_CHUNK, n_features)
        idx = np.arange(start, stop)

        n_chunk = len(idx)

        # Shape: features x samples x species
        plus = np.repeat(
            blends[None, :, :],
            n_chunk,
            axis=0
        )

        minus = plus.copy()

        for local_j, global_j in enumerate(idx):
            plus[local_j, :, global_j] += DELTA
            minus[local_j, :, global_j] -= DELTA

        plus_flat = plus.reshape(-1, n_features)
        minus_flat = minus.reshape(-1, n_features)

        p_plus = donor_probabilities(
            plus_flat, donor_idx
        ).reshape(n_chunk, N_PAIRS)

        p_minus = donor_probabilities(
            minus_flat, donor_idx
        ).reshape(n_chunk, N_PAIRS)

        effects[idx] = np.mean(
            np.abs(p_plus - p_minus),
            axis=1
        )

    return effects

# ============================================================
# 5. REPEAT ANALYSIS ACROSS RANDOM SEEDS
# ============================================================

all_results = []

for t_idx, (recipient, donor) in enumerate(TRANSITIONS):

    transition = f"{recipient} -> {donor}"

    print("\n====================================")
    print("Transition:", transition)
    print("====================================")

    for run in range(N_RUNS):

        seed = SEED + t_idx * 100000 + run

        effects = calculate_influence(
            recipient, donor, seed
        )

        # Rank 1 = highest influence
        ranks = pd.Series(effects).rank(
            ascending=False,
            method="min"
        ).to_numpy()

        for j, species in enumerate(features):

            all_results.append({
                "Transition": transition,
                "Run": run + 1,
                "Seed": seed,
                "Species": species,
                "Influence": effects[j],
                "Rank": ranks[j],
                "Top10": ranks[j] <= TOP_K
            })

        if (run + 1) % 5 == 0:
            print(f"Completed {run + 1}/{N_RUNS}")

results = pd.DataFrame(all_results)

results.to_csv(
    OUT / "Species_Influence_All_Runs.csv",
    index=False
)

# ============================================================
# 6. SPECIES-LEVEL RANK STABILITY SUMMARY
# ============================================================

summary = (
    results.groupby(["Transition", "Species"])
    .agg(
        Mean_influence=("Influence", "mean"),
        SD_influence=("Influence", "std"),
        Median_influence=("Influence", "median"),
        Influence_CI_low=(
            "Influence",
            lambda x: x.quantile(0.025)
        ),
        Influence_CI_high=(
            "Influence",
            lambda x: x.quantile(0.975)
        ),
        Mean_rank=("Rank", "mean"),
        Median_rank=("Rank", "median"),
        Rank_CI_low=(
            "Rank",
            lambda x: x.quantile(0.025)
        ),
        Rank_CI_high=(
            "Rank",
            lambda x: x.quantile(0.975)
        ),
        Top10_frequency=("Top10", "mean")
    )
    .reset_index()
)

summary = summary.sort_values(
    ["Transition", "Mean_influence"],
    ascending=[True, False]
)

summary.to_csv(
    OUT / "Species_Rank_Stability_Summary.csv",
    index=False
)

# ============================================================
# 7. RANK CONCORDANCE ACROSS RUNS
# ============================================================

concordance_rows = []

for transition in results["Transition"].unique():

    sub = results[
        results["Transition"] == transition
    ]

    rank_matrix = sub.pivot(
        index="Species",
        columns="Run",
        values="Rank"
    )

    correlations = []

    run_ids = rank_matrix.columns.tolist()

    for i in range(len(run_ids)):
        for j in range(i + 1, len(run_ids)):

            rho, _ = spearmanr(
                rank_matrix[run_ids[i]],
                rank_matrix[run_ids[j]]
            )

            correlations.append(rho)

    concordance_rows.append({
        "Transition": transition,
        "Mean_Spearman": np.nanmean(correlations),
        "Median_Spearman": np.nanmedian(correlations),
        "Min_Spearman": np.nanmin(correlations),
        "Max_Spearman": np.nanmax(correlations)
    })

concordance = pd.DataFrame(concordance_rows)

concordance.to_csv(
    OUT / "Rank_Concordance.csv",
    index=False
)

# ============================================================
# 8. TOP-10 OVERLAP ACROSS RUNS
# ============================================================

overlap_rows = []

for transition in results["Transition"].unique():

    sub = results[
        results["Transition"] == transition
    ]

    top_sets = {
        run: set(
            df.nsmallest(TOP_K, "Rank")["Species"]
        )
        for run, df in sub.groupby("Run")
    }

    overlaps = []

    run_ids = list(top_sets)

    for i in range(len(run_ids)):
        for j in range(i + 1, len(run_ids)):

            a = top_sets[run_ids[i]]
            b = top_sets[run_ids[j]]

            overlaps.append(len(a & b) / TOP_K)

    overlap_rows.append({
        "Transition": transition,
        "Mean_Top10_Overlap": np.mean(overlaps),
        "Median_Top10_Overlap": np.median(overlaps),
        "Min_Top10_Overlap": np.min(overlaps),
        "Max_Top10_Overlap": np.max(overlaps)
    })

overlap_df = pd.DataFrame(overlap_rows)

overlap_df.to_csv(
    OUT / "Top10_Overlap.csv",
    index=False
)

# ============================================================
# 9. PLOT TOP SPECIES SELECTION FREQUENCY
# ============================================================

for transition in summary["Transition"].unique():

    sub = summary[
        summary["Transition"] == transition
    ].copy()

    top = sub.nlargest(
        15, "Top10_frequency"
    ).sort_values("Top10_frequency")

    fig, ax = plt.subplots(figsize=(8, 6))

    ax.barh(
        top["Species"],
        top["Top10_frequency"] * 100
    )

    ax.set_xlabel(
        "Frequency among top 10 species (%)"
    )

    ax.set_ylabel("Species")

    ax.set_xlim(0, 100)

    ax.set_title(
        f"Species rank stability\n{transition}"
    )

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    plt.tight_layout()

    filename = transition.replace(" -> ", "_to_")

    plt.savefig(
        OUT / f"Top10_Stability_{filename}.png",
        dpi=600,
        bbox_inches="tight"
    )

    plt.show()

# ============================================================
# 10. PRINT RESULTS
# ============================================================

print("\n=== RANK CONCORDANCE ===")
print(concordance.round(4).to_string(index=False))

print("\n=== TOP-10 OVERLAP ===")
print(overlap_df.round(4).to_string(index=False))

for transition in summary["Transition"].unique():

    print("\nTop species:", transition)

    top = summary[
        summary["Transition"] == transition
    ].head(10)

    print(
        top[[
            "Species",
            "Mean_influence",
            "Median_rank",
            "Top10_frequency"
        ]].round(5).to_string(index=False)
    )

print("\nAnalysis 3 completed.")
print("Outputs saved to:", OUT.resolve())



# ============================================================
# ANALYSIS 4
# OUT-OF-FOLD CANCER PROBABILITY VALIDATION
# CRC MDT-TM FRAMEWORK
# ============================================================

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from pathlib import Path
from sklearn.pipeline import Pipeline
from sklearn.impute import SimpleImputer
from sklearn.svm import SVC
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import (
    confusion_matrix,
    balanced_accuracy_score,
    roc_auc_score
)

# ============================================================
# 1. CONFIGURATION
# ============================================================

OUT = Path("MDT_Analysis4_OOF")
OUT.mkdir(exist_ok=True)

SEED = 42
N_SPLITS = 5

LOWER_CUTOFF = 0.70
UPPER_CUTOFF = 0.80

CLASS_ORDER = ["Healthy", "Adenoma", "Cancer"]

# ============================================================
# 2. LOAD AND ALIGN DATA
# ============================================================

X = pd.read_csv(
    "TPM_clr_batch_corrected_v2.csv",
    index_col=0
)

meta = pd.read_csv(
    "CRC_metadata_aligned_bact.csv"
)

X.index = X.index.astype(str)
meta["Sample_ID"] = meta["Sample_ID"].astype(str)

assert X.index.is_unique
assert meta["Sample_ID"].is_unique

common = X.index.intersection(meta["Sample_ID"])

meta = meta.loc[
    meta["Sample_ID"].isin(common)
].copy()

X = X.loc[meta["Sample_ID"]]

assert np.array_equal(
    X.index.to_numpy(),
    meta["Sample_ID"].to_numpy()
)

X = X.apply(pd.to_numeric, errors="raise")

assert np.isfinite(X.to_numpy()).all()

y = meta["Group"].to_numpy()

print("Samples:", X.shape[0])
print("Features:", X.shape[1])
print("\nClass distribution:")
print(pd.Series(y).value_counts())

assert X.shape[0] == 995
assert set(y) == set(CLASS_ORDER)

# ============================================================
# 3. ORIGINAL SVM OBSERVER
# ============================================================

def build_observer():

    return Pipeline([
        (
            "imputer",
            SimpleImputer(strategy="median")
        ),
        (
            "svm",
            SVC(
                kernel="rbf",
                C=10,
                gamma="scale",
                probability=True,
                class_weight="balanced",
                random_state=42
            )
        )
    ])

# ============================================================
# 4. FULL-DATA PROBABILITIES
# ============================================================

full_model = build_observer()
full_model.fit(X, y)

full_probs = full_model.predict_proba(X)

full_classes = full_model.named_steps["svm"].classes_

full_prob_df = pd.DataFrame(
    full_probs,
    columns=full_classes,
    index=X.index
)

# ============================================================
# 5. OUT-OF-FOLD PROBABILITIES
# ============================================================

cv = StratifiedKFold(
    n_splits=N_SPLITS,
    shuffle=True,
    random_state=SEED
)

oof_probs = np.full(
    (len(X), len(full_classes)),
    np.nan
)

fold_assignments = np.zeros(len(X), dtype=int)

for fold, (train_idx, test_idx) in enumerate(
    cv.split(X, y), start=1
):

    print(f"Running fold {fold}/{N_SPLITS}")

    model = build_observer()

    X_train = X.iloc[train_idx]
    X_test = X.iloc[test_idx]

    y_train = y[train_idx]

    model.fit(X_train, y_train)

    fold_probs = model.predict_proba(X_test)

    fold_classes = model.named_steps["svm"].classes_

    assert np.array_equal(
        fold_classes,
        full_classes
    )

    oof_probs[test_idx] = fold_probs
    fold_assignments[test_idx] = fold

assert np.isfinite(oof_probs).all()
assert np.all(fold_assignments > 0)

oof_prob_df = pd.DataFrame(
    oof_probs,
    columns=full_classes,
    index=X.index
)

# ============================================================
# 6. PROBABILITY CATEGORY FUNCTION
# ============================================================

def assign_category(p):

    if p < LOWER_CUTOFF:
        return "Lower"

    elif p < UPPER_CUTOFF:
        return "Intermediate"

    else:
        return "Higher"

# ============================================================
# 7. COMBINE RESULTS
# ============================================================

results = meta.copy()

results["Fold"] = fold_assignments

results["P_Cancer_Full"] = (
    full_prob_df["Cancer"].to_numpy()
)

results["P_Cancer_OOF"] = (
    oof_prob_df["Cancer"].to_numpy()
)

for cls in full_classes:

    results[f"P_{cls}_OOF"] = (
        oof_prob_df[cls].to_numpy()
    )

results["Category_Full"] = (
    results["P_Cancer_Full"]
    .apply(assign_category)
)

results["Category_OOF"] = (
    results["P_Cancer_OOF"]
    .apply(assign_category)
)

results.to_csv(
    OUT / "Cancer_Probability_Full_vs_OOF.csv",
    index=False
)

# ============================================================
# 8. DIAGNOSIS x CATEGORY TABLES
# ============================================================

CATEGORY_ORDER = [
    "Lower",
    "Intermediate",
    "Higher"
]

full_table = pd.crosstab(
    results["Group"],
    results["Category_Full"]
).reindex(
    index=CLASS_ORDER,
    columns=CATEGORY_ORDER,
    fill_value=0
)

oof_table = pd.crosstab(
    results["Group"],
    results["Category_OOF"]
).reindex(
    index=CLASS_ORDER,
    columns=CATEGORY_ORDER,
    fill_value=0
)

full_table.to_csv(
    OUT / "Full_Probability_Crosstab.csv"
)

oof_table.to_csv(
    OUT / "OOF_Probability_Crosstab.csv"
)

print("\n=== FULL-DATA CROSSTAB ===")
print(full_table)

print("\n=== OUT-OF-FOLD CROSSTAB ===")
print(oof_table)

# ============================================================
# 9. ADENOMA + CANCER ONLY
# ============================================================

ac = results.loc[
    results["Group"].isin(
        ["Adenoma", "Cancer"]
    )
].copy()

ac_table = pd.crosstab(
    ac["Group"],
    ac["Category_OOF"]
).reindex(
    index=["Adenoma", "Cancer"],
    columns=CATEGORY_ORDER,
    fill_value=0
)

ac_table.to_csv(
    OUT / "OOF_Adenoma_Cancer_Crosstab.csv"
)

print("\n=== ADENOMA-CANCER OOF CATEGORIES ===")
print(ac_table)

# ============================================================
# 10. PROBABILITY SUMMARY BY DIAGNOSIS
# ============================================================

prob_summary = (
    results.groupby("Group")[
        ["P_Cancer_Full", "P_Cancer_OOF"]
    ]
    .agg(["mean", "median", "std", "min", "max"])
)

prob_summary.to_csv(
    OUT / "Cancer_Probability_Summary.csv"
)

print("\n=== PROBABILITY SUMMARY ===")
print(prob_summary.round(4))

# ============================================================
# 11. OOF CLASSIFICATION PERFORMANCE
# ============================================================

oof_labels = full_classes[
    np.argmax(oof_probs, axis=1)
]

balanced_acc = balanced_accuracy_score(
    y, oof_labels
)

print("\n=== OOF PERFORMANCE ===")
print("Balanced accuracy:", round(balanced_acc, 4))

for cls in CLASS_ORDER:

    binary_y = (y == cls).astype(int)

    auc = roc_auc_score(
        binary_y,
        oof_prob_df[cls].to_numpy()
    )

    print(f"{cls} AUROC: {auc:.4f}")

cm = confusion_matrix(
    y,
    oof_labels,
    labels=CLASS_ORDER
)

cm_df = pd.DataFrame(
    cm,
    index=CLASS_ORDER,
    columns=CLASS_ORDER
)

cm_df.to_csv(
    OUT / "OOF_Confusion_Matrix.csv"
)

print("\nOOF confusion matrix:")
print(cm_df)

# ============================================================
# 12. INTERMEDIATE BAND COUNTS
# ============================================================

intermediate_full = (
    results["Category_Full"] == "Intermediate"
).sum()

intermediate_oof = (
    results["Category_OOF"] == "Intermediate"
).sum()

intermediate_ac = (
    ac["Category_OOF"] == "Intermediate"
).sum()

print("\n=== INTERMEDIATE PROBABILITY BAND ===")
print("Full-data, all samples:", intermediate_full)
print("OOF, all samples:", intermediate_oof)
print("OOF, Adenoma + Cancer:", intermediate_ac)

# ============================================================
# 13. VISUALIZATION: OOF PROBABILITIES
# ============================================================

colors = {
    "Healthy": "#ffcc00",
    "Adenoma": "#984ea3",
    "Cancer": "#4daf4a"
}

fig, ax = plt.subplots(figsize=(8, 5))

rng = np.random.default_rng(SEED)

for i, group in enumerate(CLASS_ORDER):

    sub = results.loc[
        results["Group"] == group,
        "P_Cancer_OOF"
    ].to_numpy()

    jitter = rng.normal(
        loc=i,
        scale=0.08,
        size=len(sub)
    )

    ax.scatter(
        jitter,
        sub,
        alpha=0.35,
        s=18,
        color=colors[group]
    )

ax.axhline(
    LOWER_CUTOFF,
    color="black",
    linestyle="--",
    linewidth=1
)

ax.axhline(
    UPPER_CUTOFF,
    color="black",
    linestyle=":",
    linewidth=1
)

ax.set_xticks(range(3))
ax.set_xticklabels(CLASS_ORDER)

ax.set_ylabel("Out-of-fold P(Cancer)")
ax.set_title(
    "Cancer-class probabilities from held-out folds"
)

ax.set_ylim(-0.02, 1.02)

ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

plt.tight_layout()

plt.savefig(
    OUT / "OOF_Cancer_Probability_Distribution.png",
    dpi=600,
    bbox_inches="tight"
)

plt.show()

# ============================================================
# 14. VISUALIZATION: FULL VS OOF
# ============================================================

fig, ax = plt.subplots(figsize=(6, 6))

for group in CLASS_ORDER:

    sub = results[
        results["Group"] == group
    ]

    ax.scatter(
        sub["P_Cancer_Full"],
        sub["P_Cancer_OOF"],
        label=group,
        color=colors[group],
        alpha=0.5,
        s=20
    )

ax.plot(
    [0, 1],
    [0, 1],
    linestyle="--",
    color="black"
)

ax.axhline(
    LOWER_CUTOFF,
    linestyle=":",
    color="gray"
)

ax.axhline(
    UPPER_CUTOFF,
    linestyle=":",
    color="gray"
)

ax.set_xlabel("Full-data P(Cancer)")
ax.set_ylabel("Out-of-fold P(Cancer)")
ax.set_title("Full-data vs out-of-fold probabilities")

ax.legend(frameon=False)

ax.set_xlim(0, 1)
ax.set_ylim(0, 1)

ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

plt.tight_layout()

plt.savefig(
    OUT / "Full_vs_OOF_Cancer_Probabilities.png",
    dpi=600,
    bbox_inches="tight"
)

plt.show()

# ============================================================
# 15. FINAL SUMMARY
# ============================================================

print("\n======================================")
print("ANALYSIS 4 COMPLETED")
print("======================================")

print("Total samples:", len(results))
print("Adenoma + Cancer:", len(ac))

print("\nMain output files:")
print("- Full_Probability_Crosstab.csv")
print("- OOF_Probability_Crosstab.csv")
print("- OOF_Adenoma_Cancer_Crosstab.csv")
print("- Cancer_Probability_Summary.csv")
print("- Cancer_Probability_Full_vs_OOF.csv")
print("- OOF_Confusion_Matrix.csv")