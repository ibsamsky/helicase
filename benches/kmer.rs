use criterion::{BatchSize, Criterion, Throughput, black_box, criterion_group, criterion_main};

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

    let kmer = helicase::small::Kmer::<K>::new();
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

    let kmer = helicase::small::Kmer::<K>::new();
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

criterion_group!(benches, small, unbounded);
criterion_main!(benches);
