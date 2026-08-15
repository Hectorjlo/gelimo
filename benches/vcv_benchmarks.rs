use criterion::{black_box, criterion_group, criterion_main, Criterion};
use gelimo::data::{ComparativeData, ComparativeDataOptions};
use gelimo::math::{
    compute_vcv, solve_gls, solve_gls_linear, LinearSolverWorkspace,
};
use gelimo::opt::{grid_search_lambda_parallel, optimize_pagel_lambda};

fn bench_gls_solvers(c: &mut Criterion) {
    let options = ComparativeDataOptions {
        prune: true,
        fuzzy_taxonomy: true,
    };
    let comp_data = ComparativeData::from_files_with_options(
        "user/test/mammals_ng/tree/SpeciesTree_rooted.txt",
        "user/test/AnAge/anage_data_log_mammals_filtered.tsv",
        &options,
    )
    .unwrap();
    let y = comp_data.to_response_vector("Maximum longevity (yrs)").unwrap();
    let x = comp_data.to_design_matrix(&["Log Adult weight (g)"]).unwrap();
    let vcv = compute_vcv(&comp_data.tree);
    let mut ws = LinearSolverWorkspace::new(comp_data.tree.nodes.len(), 2);

    let mut group = c.benchmark_group("gls_solvers");

    group.bench_function("dense_cholesky_O(N3)", |b| {
        b.iter(|| {
            solve_gls(black_box(&vcv), black_box(&x), black_box(&y)).unwrap();
        })
    });

    group.bench_function("linear_ho_ane_O(N)_zero_alloc", |b| {
        b.iter(|| {
            solve_gls_linear(
                black_box(&comp_data.tree),
                black_box(&x),
                black_box(&y),
                Some(&mut ws),
            )
            .unwrap();
        })
    });

    group.finish();
}

fn bench_optimization(c: &mut Criterion) {
    let options = ComparativeDataOptions {
        prune: true,
        fuzzy_taxonomy: true,
    };
    let comp_data = ComparativeData::from_files_with_options(
        "user/test/mammals_ng/tree/SpeciesTree_rooted.txt",
        "user/test/AnAge/anage_data_log_mammals_filtered.tsv",
        &options,
    )
    .unwrap();
    let y = comp_data.to_response_vector("Maximum longevity (yrs)").unwrap();
    let x = comp_data.to_design_matrix(&["Log Adult weight (g)"]).unwrap();
    let vcv = compute_vcv(&comp_data.tree);

    let mut group = c.benchmark_group("optimization");

    group.bench_function("ml_pagel_lambda_1d_opt", |b| {
        b.iter(|| {
            optimize_pagel_lambda(
                black_box(&vcv),
                black_box(&x),
                black_box(&y),
                0.0,
                1.0,
                1e-6,
            )
            .unwrap();
        })
    });

    group.bench_function("rayon_parallel_grid_search_100_pts", |b| {
        b.iter(|| {
            grid_search_lambda_parallel(
                black_box(&vcv),
                black_box(&x),
                black_box(&y),
                100,
                0.0,
                1.0,
            );
        })
    });

    group.finish();
}

criterion_group!(benches, bench_gls_solvers, bench_optimization);
criterion_main!(benches);
