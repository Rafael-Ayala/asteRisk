# NAIF integer ID codes
# Named vector where each element is a NAIF integer ID and each name is a
# body name. Note that some NAIF integer IDs might appear multiple times
# but with different name each, to account for officially recognized synonyms
# (not including case variations)
# Official list taken from https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/req/naif_ids.html

NAIFintegerIDcodes_name2ID <- c(rep(0,3), rep(1, 2), rep(2, 2), rep(3, 5), rep(4, 2),
                                rep(5, 2), rep(6, 2), rep(7, 2), rep(8, 2), rep(9, 2),
                                10, 199, 299, 399, 301, 499, 401, 402, 599, 501:550, 553, 699,
                                601:653, 799, 701:727, 899, 801:813, 999, 901:905, -1,
                                rep(-3, 2), rep(-5, 4), rep(-6, 2), rep(-7, 2), -8,
                                rep(-12, 4), -13, rep(-18, 2), rep(-20, 2), -21,
                                rep(-23, 2), rep(-24, 2), rep(-25, 2), rep(-27, 2),
                                rep(-28, 2), rep(-29, 3), rep(-30, 3), rep(-31, 2),
                                rep(-32, 2), rep(-33, 2), rep(-37, 3), rep(-39, 2),
                                -40, rep(-41, 2), -43, rep(-44, 2), rep(-45, 2),
                                rep(-46, 2), rep(-47, 4), rep(-48, 2), -49, rep(-53, 4),
                                -55, -57, rep(-58, 2), -59, -61, rep(-62, 2), rep(-64, 2),
                                rep(-65, 2), rep(-66, 3), -67, rep(-68, 3), -70,
                                rep(-72, 2), rep(-74, 2), rep(-76, 3), rep(-77, 2),
                                -78, rep(-79, 3), -81, rep(-82, 2), -84, rep(-85, 3),
                                rep(-86, 2), -90, rep(-93, 2), rep(-94, 2), -95)