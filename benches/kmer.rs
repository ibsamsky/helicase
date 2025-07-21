use std::hint::black_box;

use criterion::{BatchSize, Criterion, Throughput, criterion_group, criterion_main};

fn bases(n: usize) -> Vec<helicase::Base> {
    let mut rng = fastrand::Rng::new();
    (0..)
        .map(|_| unsafe { helicase::Base::from_u8_unchecked(rng.u8(0..=3)) })
        .take(n)
        .collect()
}

fn small(c: &mut Criterion) {
    const K: usize = 23;
    const N: usize = 128;

    let mut group = c.benchmark_group("fixed (small) kmer");

    let mut kmer = helicase::small::Kmer::<K>::new();
    group.throughput(Throughput::Elements(N as u64));
    group.bench_function(format!("push {N} k={K}"), |b| {
        b.iter_batched(
            || bases(N),
            |bases| {
                for base in bases {
                    kmer.push(black_box(base));
                }
            },
            BatchSize::SmallInput,
        )
    });

    let mut kmer = helicase::small::Kmer::<K>::new();
    for base in std::iter::repeat_n(helicase::Base::A, K) {
        kmer.push(black_box(base));
    }

    group.throughput(Throughput::Elements(K as u64));
    group.bench_function(format!("iter {K} k={K}"), |b| {
        b.iter(|| {
            for base in kmer.bases() {
                black_box(base);
            }
        })
    });

    group.throughput(Throughput::Elements(K as u64));
    group.bench_function(format!("collect {K} k={K}"), |b| {
        b.iter_with_large_drop(|| {
            let bases: Vec<helicase::Base> = kmer.bases().collect();
            black_box(bases);
        })
    });

    group.finish();
}

#[cfg(feature = "bitvec")]
fn unbounded(c: &mut Criterion) {
    let mut group = c.benchmark_group("fixed (unbounded) kmer");

    for pow in 3..8 {
        let k = 1 << pow;
        let n = k * 2;

        let mut kmer = helicase::unbounded::Kmer::new(k);
        group.throughput(Throughput::Elements(k as u64));
        group.bench_function(format!("push {n} k={k}"), |b| {
            b.iter_batched(
                || bases(n),
                |bases| {
                    for base in bases {
                        kmer.push(black_box(base));
                    }
                },
                BatchSize::SmallInput,
            )
        });

        let mut kmer = helicase::unbounded::Kmer::new(k);
        for base in std::iter::repeat_n(helicase::Base::A, k) {
            kmer.push(black_box(base));
        }

        group.throughput(Throughput::Elements(k as u64));
        group.bench_function(format!("iter {k} k={k}"), |b| {
            b.iter(|| {
                for base in kmer.bases() {
                    black_box(base);
                }
            })
        });

        group.throughput(Throughput::Elements(k as u64));
        group.bench_function(format!("collect {k} k={k}"), |b| {
            b.iter_with_large_drop(|| {
                let bases: Vec<helicase::Base> = kmer.bases().collect();
                black_box(bases);
            })
        });
    }

    group.finish();
}

#[cfg(not(feature = "bitvec"))]
fn unbounded(_c: &mut Criterion) {}

fn sequence(c: &mut Criterion) {
    let mut group = c.benchmark_group("sequence");
    const K: usize = 23;

    for pow in [5, 10, 20] {
        let n = 1 << pow;

        group.throughput(Throughput::Elements(n as u64));
        group.bench_function(format!("push {n} bases"), |b| {
            b.iter_batched(
                || bases(n),
                |bases| {
                    let mut seq = helicase::Sequence::new();
                    for base in bases {
                        seq.push(black_box(base));
                    }
                },
                BatchSize::SmallInput,
            );
        });

        let mut seq = helicase::Sequence::new();
        for base in std::iter::repeat_n(helicase::Base::A, n) {
            seq.push(black_box(base));
        }

        let num_kmers = n - K + 1;
        group.throughput(Throughput::Elements(num_kmers as u64));
        group.bench_function(format!("iter {num_kmers} kmers k={K}"), |b| {
            b.iter(|| {
                for kmer in seq.kmers::<K>() {
                    black_box(kmer);
                }
            })
        });

        group.throughput(Throughput::Elements(num_kmers as u64));
        group.bench_function(format!("iter and mask {num_kmers} kmers k={K}"), |b| {
            b.iter(|| {
                for kmer in seq.kmers::<K>().map(|kmer| kmer.as_masked()) {
                    black_box(kmer);
                }
            })
        });
    }
}

criterion_group!(benches, small, unbounded, sequence);
criterion_main!(benches);
