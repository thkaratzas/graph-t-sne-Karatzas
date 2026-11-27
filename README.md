# Graph t-SNE με Barnes–Hut

Εργασία 1, Παράλληλα και Διανεμημένα Συστήματα Υπολογιστών, ΑΠΘ, 2026
Αθανάσιος Καρατζάς (ΑΕΜ 10514)

Υπολογισμός μίας κλίσης του t-SNE για μεγάλους αραιούς γράφους: πιθανότητες P,
έλξη κατά μήκος των ακμών (αραιό × πυκνό γινόμενο), άπωση με quadtree και
Barnes–Hut, και η κλίση 4(F_attr − F_rep). Κάθε βήμα υπάρχει σε τέσσερις εκδόσεις:
σειριακή, OpenMP, OpenCilk και Pthreads.

- Αναφορά: [`report/report-gr.pdf`](report/report-gr.pdf)
- Εκφώνηση: [`STATEMENT.md`](STATEMENT.md)

## Δομή

| Φάκελος | Περιεχόμενο |
|---|---|
| `src/` | κώδικας: `*_seq.c`, `*_omp.c`, `*_cilk.c`, `*_pt.c` και κοινά αρχεία (`tsne.h`, `p_common.c`, `quadtree.c`, `bh_common.c`, `io.c`) |
| `starter/` | ο reader των δεδομένων (`data.h`, `data.c`) από το template |
| `data/` | οι τρεις γράφοι, οι αρχικές θέσεις `-Y0` και οι τιμές αναφοράς `-ref` |
| `scripts/` | μετρήσεις και γραφήματα |
| `results/` | χρόνοι, έλεγχοι ορθότητας, στοιχεία μηχανήματος, αποτελέσματα της σειριακής έκδοσης |
| `figures/` | τα γραφήματα της αναφοράς (PDF και PNG) |
| `report/` | η αναφορά και η πρόβλεψη για το coPapersDBLP |

## Compile

Σε Linux (ή στο Codespace του template, μετά το `bash env/setup.sh`):

```bash
make
```

Φτιάχνει τα `bin/seq`, `bin/omp`, `bin/cilk`, `bin/pthreads` και `bin/compare`.

Σε macOS (Homebrew `gcc`, `libmatio`, `pkg-config` και OpenCilk 3.0) οι compilers
δηλώνονται σε ένα τοπικό `config.mk`, που δεν ανεβαίνει στο repo:

```make
CC     = /opt/homebrew/bin/gcc-16
CILKCC = $(HOME)/opencilk/bin/clang
```

Με `make info` φαίνονται οι compilers και τα flags.

## Εκτέλεση

```bash
./bin/seq data/k49_norm_10NN            # m = 16 (προεπιλογή)
./bin/seq data/k49_norm_10NN 8          # άλλη χωρητικότητα φύλλου

OMP_NUM_THREADS=4 OMP_SCHEDULE="dynamic,256" ./bin/omp data/coPapersDBLP
CILK_NWORKERS=4 ./bin/cilk data/coPapersDBLP
PT_NUM_THREADS=4 PT_SCHEDULE=dynamic PT_CHUNK=256 ./bin/pthreads data/coPapersDBLP
```

Κάθε πρόγραμμα τυπώνει τον χρόνο κάθε φάσης και το σφάλμα σε σχέση με τις τιμές
αναφοράς, και γράφει τα β, F_attr, F_rep, Z και την κλίση στο
`results/<graph>-<version>.bin`. Ο άμεσος υπολογισμός O(n²) τρέχει μόνο για γράφους
έως 50 000 κόμβους.

## Αναπαραγωγή των αποτελεσμάτων της αναφοράς

Από τη ρίζα του repo, με αυτή τη σειρά:

| Εντολή | Βγάζει | Στην αναφορά |
|---|---|---|
| `bash scripts/sweep_m.sh` | `results/sweep_m.csv` | επιλογή του m |
| `bash scripts/check.sh` | `results/check-<host>.txt` | ενότητα 4 (ορθότητα, δέντρο, συμφωνία εκδόσεων) |
| `bash scripts/measure.sh` | `results/times-<host>.csv`, `results/machine-<host>.txt` | ενότητες 5 και 6 |
| `python3 scripts/plots.py results/times-<host>.csv` | `figures/*.pdf`, `figures/*.png` | σχήματα 1–3 |

Το `sweep_m.sh` ξαναγράφει τα `results/*-seq.bin` με το τελευταίο m που δοκιμάζει,
γι' αυτό το `check.sh` τρέχει μετά από αυτό. Το `plots.py` χρειάζεται `pandas` και
`matplotlib`.

Οι μετρήσεις της αναφοράς έγιναν σε Apple M1 (4 performance + 4 efficiency πυρήνες,
8 GB, macOS 26.6.2, gcc 16.2.0 και OpenCilk 3.0, `-O2`), με τον ελάχιστο χρόνο
από 3 εκτελέσεις.

## Χρήση AI

Χρησιμοποιήθηκε το Claude (Anthropic). Για τι και τι έκανε λάθος περιγράφεται
στην ενότητα 8 της αναφοράς.
