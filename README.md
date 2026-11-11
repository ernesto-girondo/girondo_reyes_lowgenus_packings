# Extremal disc packings in surfaces of minimal genus — programs and data

*Version of 31 October 2026.* Archived at Zenodo: [doi:10.5281/zenodo.24136024](https://doi.org/10.5281/zenodo.24136024) (all versions).

*This file and `VERIFYING.md` are written in Markdown, a plain-text format that GitHub, GitLab and most editors display formatted. A formatted copy of each is included as `README.pdf` and `VERIFYING.pdf`; the `.md` files are the reference versions.*

Programs and data for the paper by E. Girondo and C. Reyes-Monsalve, *Extremal disc packings in compact hyperbolic surfaces
of minimal genus*. They determine all the extremal k-packings of all the k-extremal surfaces of genus 2
(orientable) and genus 3 (non-orientable), 423 surfaces in total, every packing being confirmed by an exact algebraic
identity.

**To check the results of the paper, read [`VERIFYING.md`](VERIFYING.md).** It explains how to install the software,
how to check any statement about a given surface, and how to recompute everything from scratch.

This file describes the contents of the folder.

## Requirements

Python 3.10 or later with `mpmath` and `cypari2` (and `numpy`, `matplotlib` for the figures, `openpyxl` for the
spreadsheet); GAP 4.12 or later.
About 2 GB of memory. All the programs are run from this folder. Installation: Section 2 of `VERIFYING.md`.

## Programs

Core
* `field.py`: the number field L_N = Q(u, cos(π/N)), u = e^R, cosh R = 1/(2 sin(π/N)), and the exact matrices of the
  reflections a, b, c.
* `geom.py`: hyperbolic geometry at 40 digits (mpmath): reduction of an isometry to a word in a, b, c, and the
  enumeration of the centres of the tessellation (admissible distances).
* `gapiface.py`: calls GAP to obtain the action of a, b, c on the cosets of a subgroup given by generators.
* `surface.py`: one surface D/K: fundamental domain and side pairings, admissible distances, candidates (intersections
  of the ED-cycles, or equi-displaced cycles, of two side pairings), search of the matrices α, exact verification, packings.
  Parameters: TOL = 1e-8, MARGIN = 0.1 (enlargement of the search disc beyond the circumradius).
* `batch.py`: runs `surface.py` on whole families, including the canonical double covers of the non-orientable
  surfaces, and computes the isometries induced by the triangle group.

Enumeration of the groups
* `enumerate.g` (GAP) and `run_enum.sh`: all torsion-free subgroups of the required index, up to conjugation.

Checks and analysis (each one writes the file indicated in `groups/`)
* `verify_packings.py` → `verify_packings.txt`: independent exact verification of all the stored packings.
* `summary_tables.py` → `surfaces_summary.csv`, `surfaces_summary.xlsx`, `double_covers.txt`: the main results in two tables that can be
  read without Python (a spreadsheet, a text editor); they are described in the next section.
* `report.py` → `report.txt`: the tables of Section 6 of the paper and consistency checks.
* `isometry_orbits.py` → `isometry_orbits.txt`: which packings are related by isometries; which classes give
  isometric surfaces.
* `identify2.py` → `surface_identifications.json`, `surface_identifications.txt`: hidden conformal and anticonformal
  automorphisms of the orientable surfaces.
* `aut_groups.py` → `automorphism_groups.txt`, `automorphism_groups.csv`: for every surface, the isomorphism type of
  its isometry group Isom(X) and, if it is orientable, of its group Aut⁺(X) of conformal automorphisms, hidden
  automorphisms included (Remark 6.9 of the paper). GAP identifies each group from its multiplication table and
  uses its own notation: C2 is the cyclic group of order 2, D12 the dihedral group of order 12 (D_6 in the paper).
  About 15 minutes.
* `genus4_covers.py` → `genus4_covers.txt`: the 144 non-orientable extremal surfaces of genus 4 of the external test
  (`nec_N18`) and their canonical double covers, 2-extremal surfaces of genus 3: which packings of the double cover
  are lifts (invariant under the covering involution), the classes of the groups K_j and the orbits of the isometry
  group (Remark 6.10 of the paper). About 10 minutes.
* `conformal_aut.py` → `conformal_aut_fuchs_N9.txt`, `conformal_aut_fuchs_N18.txt`: Aut⁺(X) for the families
  `fuchs_N9` and `fuchs_N18` only (Theorem 6.3(iv), Remarks 6.7 and 6.8), one line per surface, comparing the orders
  of its elements with those of the elements induced by the triangle group.
* `fullaut.py` → `full_aut_4_5_packings.txt`: for the surfaces with 4 or 5 packings, how each packing contributes to
  the isometry group (Theorem 6.5(i), Remark 6.6), and the fixed points of the involutions of X_6.
* `multi_extremal.py` → `multi_extremal.txt`: the 1-extremal surfaces that are also 4-extremal.
* `cross_check_9_18.py` → `cross_check_9_18.txt`: the five 4-packings of X_6 seen from the polygons of its 1-extremal tessellation (the 18-gons).
* `weierstrass.py` → `weierstrass_X6.txt`: Weierstrass points of X_6, and the monodromy of the quotient X_6/Aut⁺(X_6)
  over D/Δ(2,3,9) (Theorem 6.5(ii)).
* `margins.py` → `margins.txt`: the safety margins of the floating point decisions of the search (Section 4 of the
  paper, "Exact and floating point steps"); about 8 minutes.
* `identify_class.py`: identifies a surface given by its monodromy (the permutations of a, b, c, or the monodromy
  x, y of a dessin of type (2,3,N)) and finds it in the catalogue, whatever the numbering of the sheets; e.g.
  `python3 identify_class.py 18 xy "<x>" "<y>"`. Plain Python, no PARI or GAP. Section 4.7 of `VERIFYING.md`.
* `check_surface.py`: runs the search again for one surface and prints its steps; with `--details`, every admissible
  distance, the two side pairings of Step 3 and every intersection point of their ED-cycles with its fate
  (discarded, centre of a disc of a packing, or shown not to be a centre). Section 4.3 of `VERIFYING.md`.
* `autfix.py`: explicit automorphisms, their orders and fixed points, for any surface (prints only).
* `draw_surface.py` → folder `surface_figures/`: a picture of any surface (fundamental domain, side pairings, centres
  of all its extremal packings; fixed points and mirrors of all its isometries), e.g.
  `python3 draw_surface.py fuchs_N9 21`. Uses `isometries.py`. With `--edcycles`, a picture of Step 3 of the search
  (the ED-cycles of the two side pairings, all their intersection points in the search disc and the candidates, as
  Figure 2 of the paper does for X_6), with a text file listing every point; the picture can be crowded when many
  ED-cycles overlap (typically for N = 7, 8), the text file always has the complete information.
* `dirichlet.py`: Dirichlet domains centred at the disc centres, with their side pairings and checks (vertex
  cycles, area); `python3 draw_surface.py fuchs_N9 21 --dirichlet-all` draws them.
* `draw_all.py` → folder `surface_figures/all/`: the pictures (packings, isometries, Dirichlet domains, ED-cycles of
  Step 3 with their text files) of all the
  373 surfaces, with a table and a consistency check against the data files (`checks.txt`).
* `figures.py` → folder `figures/`: the figures of the paper (Figures 4–7 use `draw_surface.py` and `dirichlet.py`).
* `regenerate_all.sh`, `compare_outputs.py`: recompute everything and compare with the distributed data.
* `wordfmt.py`: writes words in a, b, c with powers, `(ca)^3cb` for `cacacacb` (used in the tables and reports).
* `outfile.py`: small helper that copies the printed output of a program to its file.

## Data (`groups/`)

* `<family>.jsonl`: the groups computed by GAP, one per line (class #n is line n): the permutations of a, b, c on the
  cosets of K (point 1 = K) and generators of K. Families: `nec_N12`, `nec_N9`, `nec_N8`, `nec_N7` (non-orientable
  surfaces of genus 3, subgroups of index 24, 36, 48, 84 of Δ±(2,3,N)) and `fuchs_N18`, `fuchs_N12`, `fuchs_N10`,
  `fuchs_N9`, `fuchs_N8`, `fuchs_N7` (orientable surfaces of genus 2, subgroups of index 18, 24, 30, 36, 48, 84 of
  Δ(2,3,N)). `nec_N18` (index 36 in Δ±(2,3,18)) is an external test: the non-orientable extremal surfaces of genus 4.
* `<family>_results.json`: one record per surface: k, N, index, orientability, side pairings (words in a, b, c),
  number of candidates, and for each extremal packing the centre of one disc in the disc model, its distance to O,
  the exact matrix α (a PARI matrix over Q[y]/(P_N)) and the words w_s of the identity s·α = ±α·w_s; the isometries induced by
  the triangle group (coset, parity, order, fixed points); the canonical form of the class (field `key`); the number
  of matrices α that could not be verified exactly (field `n_uncertified_alphas`, always 0) and the confirmation of the rejected candidates. For the non-orientable
  surfaces, the field `double_cover` contains the same data for the canonical double cover K⁺ = K ∩ Δ.
* `surfaces_summary.xlsx` and `surfaces_summary.csv` (the same table; the `.xlsx` opens directly in Excel, Numbers or
  LibreOffice, while some spreadsheets in Spanish or other languages put the whole `.csv` line in the first column,
  because they expect `;` and not `,` between the columns): one row for each of the 423 surfaces (the 373 classes and the 50 canonical double covers,
  named `nec_N9#2+` for the double cover of `nec_N9#2`). Columns: `surface`; `N`; `k`; `orientable`;
  `double_cover_of` (only for the double covers); `index_in_Delta_pm` (index of K in Δ±(2,3,N)); `n_packings`
  (number of extremal k-packings); `induced_isometries` (order of the group of isometries of the surface induced by
  Δ±(2,3,N)); `induced_anticonformal` (whether that group contains an orientation reversing isometry);
  `centre_packing_2`, …, `centre_packing_5` (for every packing other than the one of the uniformization, the
  centre of one of its discs, in the disc model, to 8 decimals); `class_id` (a short fingerprint of the canonical
  form `key`: it does not depend on the order in which GAP lists the classes, and it is the same for two classes
  conjugate in Δ±(2,3,N), i.e. a class and its mirror image, and for a double cover and its orientable class);
  `generators_of_K` (generators of the surface group K itself, as words in a, b, c, as returned by GAP; GAP
  writes inverses such as `c^-1` because it treats a, b, c as free generators, but a, b, c are reflections, so
  `c^-1` = `c`; this column is empty for the double covers); `side_pairings` (the side pairings of
  the fundamental domain, as words in a, b, c, written with powers of repeated blocks: `(ca)^3cb(ac)^2` is
  `cacacacbacac`; the words are read from left to right).
* `double_covers.txt`: for each non-orientable surface, the class of the orientable surfaces that is its canonical
  double cover, with the number of packings and the order of the induced isometry group of both (Appendix A of the
  paper), and, for each N, the orientable classes with an induced anticonformal isometry that are not double covers.
* `<family>_run.log`, `enum.log`: logs of the computation.
* The output files of the analysis programs listed above.

## Conventions

* a, b, c are the reflections in the sides of the triangle with angles π/N (at the centre O of the polygon: i in the
  upper half plane, 0 in the disc), π/2 (midpoint of an edge) and π/3 (vertex): a in the line through O and the edge
  midpoint, b in the edge, c in the line through O and the vertex. a² = b² = c² = (ab)² = (bc)³ = (ca)^N = 1.
  Matrices act on the upper half plane, the orientation reversing ones on the complex conjugate. A word is read
  from left to right, both as a product of matrices and in the action on cosets.
* Coset actions are right actions on {1, …, d} (GAP's `FactorCosetAction`), point 1 being the coset K; in the Python
  code points are numbered from 0.
* Disc model: z ↦ i(z − i)/(z + i); the first edge midpoint of the central polygon lies on the negative imaginary
  axis. Angles in the paper are measured in this model, from this position of the central polygon.

## Validation

* The surface K_(9,−,2) worked out in detail in Section 7.1 of C. Reyes-Monsalve's thesis is reproduced: the permutations of
  a, b, c printed there define exactly the class `nec_N9#2` (with a different numbering of the cosets), and the
  program finds the same first 15 admissible distances (d_0, …, d_14), a fundamental domain with 8 side pairings as
  there, the same hidden centre (≈ 0.4475 + 0.0165i) and two packings, also two for the double cover. The
  identification can be repeated with `identify_class.py` (Example 1 of Section 4.7 of `VERIFYING.md`).
* All the class counts agree with Table 2 of the paper.
* The 50 non-orientable surfaces agree with the classification in C. Reyes-Monsalve's thesis.
* The eight 1-extremal surfaces of genus 2 agree in every detail with Girondo–González-Diez, Israel J. Math. 132
  (2002).
* The 144 non-orientable extremal surfaces of genus 4 reproduce Nakamura, Conform. Geom. Dyn. 13 (2009): 22 with
  two extremal discs and 122 with one, the positions of the second disc and the automorphism groups of his
  Tables 1–2. This test revealed a bug, since corrected (a degenerate ED-cycle, when an admissible distance
  equals the translation length of a side pairing); the results of the paper did not change.
* Every candidate rejected in the last step of the search is moved by some product of two side pairings (other
  than the identity) a distance that is not admissible, so it is not a centre; for 10 of the 70 such candidates
  (on `fuchs_N7#70` and `#71`) this distance is beyond the list of admissible distances of the search, and it is
  decided with a longer complete list (Lemma 3.1). Every matrix α found by the search was verified exactly.
* Every floating point decision of the search is taken with a wide margin: exact equalities are reproduced with an
  error below 2·10⁻¹⁹, rejected points miss every admissible distance by at least 6·10⁻⁶, and the confirmations of the rejected candidates
  by at least 2.7·10⁻³ (`margins.txt`).
* The computation is deterministic; the whole folder has been regenerated from a clean copy with identical output
  (Section 6 of `VERIFYING.md`).

## License

The programs are distributed under the MIT License and the data, pictures and documentation under the Creative
Commons Attribution 4.0 International License (CC BY 4.0); see the file `LICENSE`. The file `CITATION.cff` says how
to cite this collection.

## Credits

The method and earlier versions of the programs are due to E. Girondo and C. Reyes-Monsalve (Chapter 7 of C. Reyes-Monsalve's
thesis). Claude (Anthropic), an AI model, under their direction, organized, cleaned and restructured those programs
into the present ones, ran the computations, organized the data, wrote, at their request and to give more complete
information, the auxiliary programs (pictures, Dirichlet domains, automorphism groups, identification of a surface
from its monodromy), and wrote this documentation. The
authors take responsibility for all of it (see the section "Use of AI tools" of the paper).
