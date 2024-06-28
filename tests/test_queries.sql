
-- Meta data queries
SELECT count(DISTINCT SCP_acc) FROM scp_data;


-- Genrate the protien set tetramer counts

SELECT genome_id, length(tetramers), 0 as source_table from `PF00119.20_genomes` UNION ALL
SELECT genome_id, length(tetramers), 1 as source_table from `PF00121.18_genomes` UNION ALL 
SELECT genome_id, length(tetramers), 2 as source_table from `PF00162.19_genomes` UNION ALL 
SELECT genome_id, length(tetramers), 3 as source_table from `PF00164.25_genomes` ;



-- Generate the tetramer genome corrrespondence for generating genome pairs

SELECT tetramer, genomes, 0 as source_table FROM `PF00119.20_tetras` WHERE tetramer BETWEEN 2000 and 3000  UNION ALL
 SELECT tetramer, genomes, 1 as source_table FROM `PF00121.18_tetras` WHERE tetramer BETWEEN 2000 and 3000  UNION ALL
 SELECT tetramer, genomes, 2 as source_table FROM `PF00162.19_tetras` WHERE tetramer BETWEEN 2000 and 3000  UNION ALL
 SELECT tetramer, genomes, 3 as source_table FROM `PF00164.25_tetras` WHERE tetramer BETWEEN 2000 and 3000  UNION ALL
 SELECT tetramer, genomes, 4 as source_table FROM `PF00177.21_tetras` WHERE tetramer BETWEEN 2000 and 3000  UNION ALL
 SELECT tetramer, genomes, 5 as source_table FROM `PF00181.23_tetras` WHERE tetramer BETWEEN 2000 and 3000  UNION ALL
 SELECT tetramer, genomes, 6 as source_table FROM `PF00189.20_tetras` WHERE tetramer BETWEEN 2000 and 3000  UNION ALL
 SELECT tetramer, genomes, 7 as source_table FROM `PF00203.21_tetras` WHERE tetramer BETWEEN 2000 and 3000  UNION ALL
 SELECT tetramer, genomes, 8 as source_table FROM `PF00213.18_tetras` WHERE tetramer BETWEEN 2000 and 3000  UNION ALL
 SELECT tetramer, genomes, 9 as source_table FROM `PF00231.19_tetras` WHERE tetramer BETWEEN 2000 and 3000  UNION ALL
 SELECT tetramer, genomes, 10 as source_table FROM `PF00237.19_tetras` WHERE tetramer BETWEEN 2000 and 3000  UNION ALL
 SELECT tetramer, genomes, 11 as source_table FROM `PF00238.19_tetras` WHERE tetramer BETWEEN 2000 and 3000  UNION ALL
 SELECT tetramer, genomes, 12 as source_table FROM `PF00252.18_tetras` WHERE tetramer BETWEEN 2000 and 3000  UNION ALL
 SELECT tetramer, genomes, 13 as source_table FROM `PF00276.20_tetras` WHERE tetramer BETWEEN 2000 and 3000  UNION ALL
 SELECT tetramer, genomes, 14 as source_table FROM `PF00281.19_tetras` WHERE tetramer BETWEEN 2000 and 3000  UNION ALL
 SELECT tetramer, genomes, 15 as source_table FROM `PF00297.22_tetras` WHERE tetramer BETWEEN 2000 and 3000  UNION ALL
 SELECT tetramer, genomes, 16 as source_table FROM `PF00312.22_tetras` WHERE tetramer BETWEEN 2000 and 3000  UNION ALL
 SELECT tetramer, genomes, 17 as source_table FROM `PF00318.20_tetras` WHERE tetramer BETWEEN 2000 and 3000  UNION ALL
 SELECT tetramer, genomes, 18 as source_table FROM `PF00334.19_tetras` WHERE tetramer BETWEEN 2000 and 3000  UNION ALL
 SELECT tetramer, genomes, 19 as source_table FROM `PF00338.22_tetras` WHERE tetramer BETWEEN 2000 and 3000  UNION ALL
 SELECT tetramer, genomes, 20 as source_table FROM `PF00344.20_tetras` WHERE tetramer BETWEEN 2000 and 3000  UNION ALL
 SELECT tetramer, genomes, 21 as source_table FROM `PF00347.23_tetras` WHERE tetramer BETWEEN 2000 and 3000  UNION ALL
 SELECT tetramer, genomes, 22 as source_table FROM `PF00366.20_tetras` WHERE tetramer BETWEEN 2000 and 3000  UNION ALL
 SELECT tetramer, genomes, 23 as source_table FROM `PF00380.19_tetras` WHERE tetramer BETWEEN 2000 and 3000  UNION ALL
 SELECT tetramer, genomes, 24 as source_table FROM `PF00406.22_tetras` WHERE tetramer BETWEEN 2000 and 3000  UNION ALL
 SELECT tetramer, genomes, 25 as source_table FROM `PF00410.19_tetras` WHERE tetramer BETWEEN 2000 and 3000  UNION ALL
 SELECT tetramer, genomes, 26 as source_table FROM `PF00411.19_tetras` WHERE tetramer BETWEEN 2000 and 3000  UNION ALL
 SELECT tetramer, genomes, 27 as source_table FROM `PF00416.22_tetras` WHERE tetramer BETWEEN 2000 and 3000  UNION ALL
 SELECT tetramer, genomes, 28 as source_table FROM `PF00453.18_tetras` WHERE tetramer BETWEEN 2000 and 3000  UNION ALL
 SELECT tetramer, genomes, 29 as source_table FROM `PF00572.18_tetras` WHERE tetramer BETWEEN 2000 and 3000  UNION ALL
 SELECT tetramer, genomes, 30 as source_table FROM `PF00573.22_tetras` WHERE tetramer BETWEEN 2000 and 3000  UNION ALL
 SELECT tetramer, genomes, 31 as source_table FROM `PF00584.20_tetras` WHERE tetramer BETWEEN 2000 and 3000  UNION ALL
 SELECT tetramer, genomes, 32 as source_table FROM `PF00687.21_tetras` WHERE tetramer BETWEEN 2000 and 3000  UNION ALL
 SELECT tetramer, genomes, 33 as source_table FROM `PF00709.21_tetras` WHERE tetramer BETWEEN 2000 and 3000  UNION ALL
 SELECT tetramer, genomes, 34 as source_table FROM `PF00749.21_tetras` WHERE tetramer BETWEEN 2000 and 3000  UNION ALL
 SELECT tetramer, genomes, 35 as source_table FROM `PF00750.19_tetras` WHERE tetramer BETWEEN 2000 and 3000  UNION ALL
 SELECT tetramer, genomes, 36 as source_table FROM `PF00825.18_tetras` WHERE tetramer BETWEEN 2000 and 3000  UNION ALL
 SELECT tetramer, genomes, 37 as source_table FROM `PF00828.19_tetras` WHERE tetramer BETWEEN 2000 and 3000  UNION ALL
 SELECT tetramer, genomes, 38 as source_table FROM `PF00829.21_tetras` WHERE tetramer BETWEEN 2000 and 3000  UNION ALL
 SELECT tetramer, genomes, 39 as source_table FROM `PF00830.19_tetras` WHERE tetramer BETWEEN 2000 and 3000  UNION ALL
 SELECT tetramer, genomes, 40 as source_table FROM `PF00831.23_tetras` WHERE tetramer BETWEEN 2000 and 3000  UNION ALL
 SELECT tetramer, genomes, 41 as source_table FROM `PF00861.22_tetras` WHERE tetramer BETWEEN 2000 and 3000  UNION ALL
 SELECT tetramer, genomes, 42 as source_table FROM `PF00886.19_tetras` WHERE tetramer BETWEEN 2000 and 3000  UNION ALL
 SELECT tetramer, genomes, 43 as source_table FROM `PF00889.19_tetras` WHERE tetramer BETWEEN 2000 and 3000  UNION ALL
 SELECT tetramer, genomes, 44 as source_table FROM `PF01016.19_tetras` WHERE tetramer BETWEEN 2000 and 3000  UNION ALL
 SELECT tetramer, genomes, 45 as source_table FROM `PF01025.19_tetras` WHERE tetramer BETWEEN 2000 and 3000  UNION ALL
 SELECT tetramer, genomes, 46 as source_table FROM `PF01142.18_tetras` WHERE tetramer BETWEEN 2000 and 3000  UNION ALL
 SELECT tetramer, genomes, 47 as source_table FROM `PF01176.19_tetras` WHERE tetramer BETWEEN 2000 and 3000  UNION ALL
 SELECT tetramer, genomes, 48 as source_table FROM `PF01192.22_tetras` WHERE tetramer BETWEEN 2000 and 3000  UNION ALL
 SELECT tetramer, genomes, 49 as source_table FROM `PF01193.24_tetras` WHERE tetramer BETWEEN 2000 and 3000  UNION ALL
 SELECT tetramer, genomes, 50 as source_table FROM `PF01195.19_tetras` WHERE tetramer BETWEEN 2000 and 3000  UNION ALL
 SELECT tetramer, genomes, 51 as source_table FROM `PF01196.19_tetras` WHERE tetramer BETWEEN 2000 and 3000  UNION ALL
 SELECT tetramer, genomes, 52 as source_table FROM `PF01245.20_tetras` WHERE tetramer BETWEEN 2000 and 3000  UNION ALL
 SELECT tetramer, genomes, 53 as source_table FROM `PF01250.17_tetras` WHERE tetramer BETWEEN 2000 and 3000  UNION ALL
 SELECT tetramer, genomes, 54 as source_table FROM `PF01264.21_tetras` WHERE tetramer BETWEEN 2000 and 3000  UNION ALL
 SELECT tetramer, genomes, 55 as source_table FROM `PF01351.18_tetras` WHERE tetramer BETWEEN 2000 and 3000  UNION ALL
 SELECT tetramer, genomes, 56 as source_table FROM `PF01632.19_tetras` WHERE tetramer BETWEEN 2000 and 3000  UNION ALL
 SELECT tetramer, genomes, 57 as source_table FROM `PF01649.18_tetras` WHERE tetramer BETWEEN 2000 and 3000  UNION ALL
 SELECT tetramer, genomes, 58 as source_table FROM `PF01668.18_tetras` WHERE tetramer BETWEEN 2000 and 3000  UNION ALL
 SELECT tetramer, genomes, 59 as source_table FROM `PF01715.17_tetras` WHERE tetramer BETWEEN 2000 and 3000  UNION ALL
 SELECT tetramer, genomes, 60 as source_table FROM `PF01725.16_tetras` WHERE tetramer BETWEEN 2000 and 3000  UNION ALL
 SELECT tetramer, genomes, 61 as source_table FROM `PF01746.21_tetras` WHERE tetramer BETWEEN 2000 and 3000  UNION ALL
 SELECT tetramer, genomes, 62 as source_table FROM `PF01765.19_tetras` WHERE tetramer BETWEEN 2000 and 3000  UNION ALL
 SELECT tetramer, genomes, 63 as source_table FROM `PF01783.23_tetras` WHERE tetramer BETWEEN 2000 and 3000  UNION ALL
 SELECT tetramer, genomes, 64 as source_table FROM `PF01808.18_tetras` WHERE tetramer BETWEEN 2000 and 3000  UNION ALL
 SELECT tetramer, genomes, 65 as source_table FROM `PF02033.18_tetras` WHERE tetramer BETWEEN 2000 and 3000  UNION ALL
 SELECT tetramer, genomes, 66 as source_table FROM `PF02130.17_tetras` WHERE tetramer BETWEEN 2000 and 3000  UNION ALL
 SELECT tetramer, genomes, 67 as source_table FROM `PF02367.17_tetras` WHERE tetramer BETWEEN 2000 and 3000  UNION ALL
 SELECT tetramer, genomes, 68 as source_table FROM `PF02410.15_tetras` WHERE tetramer BETWEEN 2000 and 3000  UNION ALL
 SELECT tetramer, genomes, 69 as source_table FROM `PF02565.15_tetras` WHERE tetramer BETWEEN 2000 and 3000  UNION ALL
 SELECT tetramer, genomes, 70 as source_table FROM `PF02601.15_tetras` WHERE tetramer BETWEEN 2000 and 3000  UNION ALL
 SELECT tetramer, genomes, 71 as source_table FROM `PF02699.15_tetras` WHERE tetramer BETWEEN 2000 and 3000  UNION ALL
 SELECT tetramer, genomes, 72 as source_table FROM `PF03652.15_tetras` WHERE tetramer BETWEEN 2000 and 3000  UNION ALL
 SELECT tetramer, genomes, 73 as source_table FROM `PF03840.14_tetras` WHERE tetramer BETWEEN 2000 and 3000  UNION ALL
 SELECT tetramer, genomes, 74 as source_table FROM `PF03948.14_tetras` WHERE tetramer BETWEEN 2000 and 3000  UNION ALL
 SELECT tetramer, genomes, 75 as source_table FROM `PF05221.17_tetras` WHERE tetramer BETWEEN 2000 and 3000  UNION ALL
 SELECT tetramer, genomes, 76 as source_table FROM `PF06026.14_tetras` WHERE tetramer BETWEEN 2000 and 3000  UNION ALL
 SELECT tetramer, genomes, 77 as source_table FROM `PF13393.6_tetras` WHERE tetramer BETWEEN 2000 and 3000  UNION ALL
 SELECT tetramer, genomes, 78 as source_table FROM `PF17136.4_tetras` WHERE tetramer BETWEEN 2000 and 3000  UNION ALL
 SELECT tetramer, genomes, 79 as source_table FROM `PF01139.17_tetras` WHERE tetramer BETWEEN 2000 and 3000  ORDER BY tetramer, source_table;

-- Get the coutns for above queries

SELECT tetramer, length(genomes), 0 as source_table FROM `PF00119.20_tetras` WHERE tetramer BETWEEN 2000 and 3000  UNION ALL
 SELECT tetramer, length(genomes), 1 as source_table FROM `PF00121.18_tetras` WHERE tetramer BETWEEN 2000 and 3000  UNION ALL
 SELECT tetramer, length(genomes), 2 as source_table FROM `PF00162.19_tetras` WHERE tetramer BETWEEN 2000 and 3000  UNION ALL
 SELECT tetramer, length(genomes), 3 as source_table FROM `PF00164.25_tetras` WHERE tetramer BETWEEN 2000 and 3000  UNION ALL
 SELECT tetramer, length(genomes), 4 as source_table FROM `PF00177.21_tetras` WHERE tetramer BETWEEN 2000 and 3000  UNION ALL
 SELECT tetramer, length(genomes), 5 as source_table FROM `PF00181.23_tetras` WHERE tetramer BETWEEN 2000 and 3000  UNION ALL
 SELECT tetramer, length(genomes), 6 as source_table FROM `PF00189.20_tetras` WHERE tetramer BETWEEN 2000 and 3000  UNION ALL
 SELECT tetramer, length(genomes), 7 as source_table FROM `PF00203.21_tetras` WHERE tetramer BETWEEN 2000 and 3000  UNION ALL
 SELECT tetramer, length(genomes), 8 as source_table FROM `PF00213.18_tetras` WHERE tetramer BETWEEN 2000 and 3000  UNION ALL
 SELECT tetramer, length(genomes), 9 as source_table FROM `PF00231.19_tetras` WHERE tetramer BETWEEN 2000 and 3000  UNION ALL
 SELECT tetramer, length(genomes), 10 as source_table FROM `PF00237.19_tetras` WHERE tetramer BETWEEN 2000 and 3000  UNION ALL
 SELECT tetramer, length(genomes), 11 as source_table FROM `PF00238.19_tetras` WHERE tetramer BETWEEN 2000 and 3000  UNION ALL
 SELECT tetramer, length(genomes), 12 as source_table FROM `PF00252.18_tetras` WHERE tetramer BETWEEN 2000 and 3000  UNION ALL
 SELECT tetramer, length(genomes), 13 as source_table FROM `PF00276.20_tetras` WHERE tetramer BETWEEN 2000 and 3000  UNION ALL
 SELECT tetramer, length(genomes), 14 as source_table FROM `PF00281.19_tetras` WHERE tetramer BETWEEN 2000 and 3000  UNION ALL
 SELECT tetramer, length(genomes), 15 as source_table FROM `PF00297.22_tetras` WHERE tetramer BETWEEN 2000 and 3000  UNION ALL
 SELECT tetramer, length(genomes), 16 as source_table FROM `PF00312.22_tetras` WHERE tetramer BETWEEN 2000 and 3000  UNION ALL
 SELECT tetramer, length(genomes), 17 as source_table FROM `PF00318.20_tetras` WHERE tetramer BETWEEN 2000 and 3000  UNION ALL
 SELECT tetramer, length(genomes), 18 as source_table FROM `PF00334.19_tetras` WHERE tetramer BETWEEN 2000 and 3000  UNION ALL
 SELECT tetramer, length(genomes), 19 as source_table FROM `PF00338.22_tetras` WHERE tetramer BETWEEN 2000 and 3000  UNION ALL
 SELECT tetramer, length(genomes), 20 as source_table FROM `PF00344.20_tetras` WHERE tetramer BETWEEN 2000 and 3000  UNION ALL
 SELECT tetramer, length(genomes), 21 as source_table FROM `PF00347.23_tetras` WHERE tetramer BETWEEN 2000 and 3000  UNION ALL
 SELECT tetramer, length(genomes), 22 as source_table FROM `PF00366.20_tetras` WHERE tetramer BETWEEN 2000 and 3000  UNION ALL
 SELECT tetramer, length(genomes), 23 as source_table FROM `PF00380.19_tetras` WHERE tetramer BETWEEN 2000 and 3000  UNION ALL
 SELECT tetramer, length(genomes), 24 as source_table FROM `PF00406.22_tetras` WHERE tetramer BETWEEN 2000 and 3000  UNION ALL
 SELECT tetramer, length(genomes), 25 as source_table FROM `PF00410.19_tetras` WHERE tetramer BETWEEN 2000 and 3000  UNION ALL
 SELECT tetramer, length(genomes), 26 as source_table FROM `PF00411.19_tetras` WHERE tetramer BETWEEN 2000 and 3000  UNION ALL
 SELECT tetramer, length(genomes), 27 as source_table FROM `PF00416.22_tetras` WHERE tetramer BETWEEN 2000 and 3000  UNION ALL
 SELECT tetramer, length(genomes), 28 as source_table FROM `PF00453.18_tetras` WHERE tetramer BETWEEN 2000 and 3000  UNION ALL
 SELECT tetramer, length(genomes), 29 as source_table FROM `PF00572.18_tetras` WHERE tetramer BETWEEN 2000 and 3000  UNION ALL
 SELECT tetramer, length(genomes), 30 as source_table FROM `PF00573.22_tetras` WHERE tetramer BETWEEN 2000 and 3000  UNION ALL
 SELECT tetramer, length(genomes), 31 as source_table FROM `PF00584.20_tetras` WHERE tetramer BETWEEN 2000 and 3000  UNION ALL
 SELECT tetramer, length(genomes), 32 as source_table FROM `PF00687.21_tetras` WHERE tetramer BETWEEN 2000 and 3000  UNION ALL
 SELECT tetramer, length(genomes), 33 as source_table FROM `PF00709.21_tetras` WHERE tetramer BETWEEN 2000 and 3000  UNION ALL
 SELECT tetramer, length(genomes), 34 as source_table FROM `PF00749.21_tetras` WHERE tetramer BETWEEN 2000 and 3000  UNION ALL
 SELECT tetramer, length(genomes), 35 as source_table FROM `PF00750.19_tetras` WHERE tetramer BETWEEN 2000 and 3000  UNION ALL
 SELECT tetramer, length(genomes), 36 as source_table FROM `PF00825.18_tetras` WHERE tetramer BETWEEN 2000 and 3000  UNION ALL
 SELECT tetramer, length(genomes), 37 as source_table FROM `PF00828.19_tetras` WHERE tetramer BETWEEN 2000 and 3000  UNION ALL
 SELECT tetramer, length(genomes), 38 as source_table FROM `PF00829.21_tetras` WHERE tetramer BETWEEN 2000 and 3000  UNION ALL
 SELECT tetramer, length(genomes), 39 as source_table FROM `PF00830.19_tetras` WHERE tetramer BETWEEN 2000 and 3000  UNION ALL
 SELECT tetramer, length(genomes), 40 as source_table FROM `PF00831.23_tetras` WHERE tetramer BETWEEN 2000 and 3000  UNION ALL
 SELECT tetramer, length(genomes), 41 as source_table FROM `PF00861.22_tetras` WHERE tetramer BETWEEN 2000 and 3000  UNION ALL
 SELECT tetramer, length(genomes), 42 as source_table FROM `PF00886.19_tetras` WHERE tetramer BETWEEN 2000 and 3000  UNION ALL
 SELECT tetramer, length(genomes), 43 as source_table FROM `PF00889.19_tetras` WHERE tetramer BETWEEN 2000 and 3000  UNION ALL
 SELECT tetramer, length(genomes), 44 as source_table FROM `PF01016.19_tetras` WHERE tetramer BETWEEN 2000 and 3000  UNION ALL
 SELECT tetramer, length(genomes), 45 as source_table FROM `PF01025.19_tetras` WHERE tetramer BETWEEN 2000 and 3000  UNION ALL
 SELECT tetramer, length(genomes), 46 as source_table FROM `PF01142.18_tetras` WHERE tetramer BETWEEN 2000 and 3000  UNION ALL
 SELECT tetramer, length(genomes), 47 as source_table FROM `PF01176.19_tetras` WHERE tetramer BETWEEN 2000 and 3000  UNION ALL
 SELECT tetramer, length(genomes), 48 as source_table FROM `PF01192.22_tetras` WHERE tetramer BETWEEN 2000 and 3000  UNION ALL
 SELECT tetramer, length(genomes), 49 as source_table FROM `PF01193.24_tetras` WHERE tetramer BETWEEN 2000 and 3000  UNION ALL
 SELECT tetramer, length(genomes), 50 as source_table FROM `PF01195.19_tetras` WHERE tetramer BETWEEN 2000 and 3000  UNION ALL
 SELECT tetramer, length(genomes), 51 as source_table FROM `PF01196.19_tetras` WHERE tetramer BETWEEN 2000 and 3000  UNION ALL
 SELECT tetramer, length(genomes), 52 as source_table FROM `PF01245.20_tetras` WHERE tetramer BETWEEN 2000 and 3000  UNION ALL
 SELECT tetramer, length(genomes), 53 as source_table FROM `PF01250.17_tetras` WHERE tetramer BETWEEN 2000 and 3000  UNION ALL
 SELECT tetramer, length(genomes), 54 as source_table FROM `PF01264.21_tetras` WHERE tetramer BETWEEN 2000 and 3000  UNION ALL
 SELECT tetramer, length(genomes), 55 as source_table FROM `PF01351.18_tetras` WHERE tetramer BETWEEN 2000 and 3000  UNION ALL
 SELECT tetramer, length(genomes), 56 as source_table FROM `PF01632.19_tetras` WHERE tetramer BETWEEN 2000 and 3000  UNION ALL
 SELECT tetramer, length(genomes), 57 as source_table FROM `PF01649.18_tetras` WHERE tetramer BETWEEN 2000 and 3000  UNION ALL
 SELECT tetramer, length(genomes), 58 as source_table FROM `PF01668.18_tetras` WHERE tetramer BETWEEN 2000 and 3000  UNION ALL
 SELECT tetramer, length(genomes), 59 as source_table FROM `PF01715.17_tetras` WHERE tetramer BETWEEN 2000 and 3000  UNION ALL
 SELECT tetramer, length(genomes), 60 as source_table FROM `PF01725.16_tetras` WHERE tetramer BETWEEN 2000 and 3000  UNION ALL
 SELECT tetramer, length(genomes), 61 as source_table FROM `PF01746.21_tetras` WHERE tetramer BETWEEN 2000 and 3000  UNION ALL
 SELECT tetramer, length(genomes), 62 as source_table FROM `PF01765.19_tetras` WHERE tetramer BETWEEN 2000 and 3000  UNION ALL
 SELECT tetramer, length(genomes), 63 as source_table FROM `PF01783.23_tetras` WHERE tetramer BETWEEN 2000 and 3000  UNION ALL
 SELECT tetramer, length(genomes), 64 as source_table FROM `PF01808.18_tetras` WHERE tetramer BETWEEN 2000 and 3000  UNION ALL
 SELECT tetramer, length(genomes), 65 as source_table FROM `PF02033.18_tetras` WHERE tetramer BETWEEN 2000 and 3000  UNION ALL
 SELECT tetramer, length(genomes), 66 as source_table FROM `PF02130.17_tetras` WHERE tetramer BETWEEN 2000 and 3000  UNION ALL
 SELECT tetramer, length(genomes), 67 as source_table FROM `PF02367.17_tetras` WHERE tetramer BETWEEN 2000 and 3000  UNION ALL
 SELECT tetramer, length(genomes), 68 as source_table FROM `PF02410.15_tetras` WHERE tetramer BETWEEN 2000 and 3000  UNION ALL
 SELECT tetramer, length(genomes), 69 as source_table FROM `PF02565.15_tetras` WHERE tetramer BETWEEN 2000 and 3000  UNION ALL
 SELECT tetramer, length(genomes), 70 as source_table FROM `PF02601.15_tetras` WHERE tetramer BETWEEN 2000 and 3000  UNION ALL
 SELECT tetramer, length(genomes), 71 as source_table FROM `PF02699.15_tetras` WHERE tetramer BETWEEN 2000 and 3000  UNION ALL
 SELECT tetramer, length(genomes), 72 as source_table FROM `PF03652.15_tetras` WHERE tetramer BETWEEN 2000 and 3000  UNION ALL
 SELECT tetramer, length(genomes), 73 as source_table FROM `PF03840.14_tetras` WHERE tetramer BETWEEN 2000 and 3000  UNION ALL
 SELECT tetramer, length(genomes), 74 as source_table FROM `PF03948.14_tetras` WHERE tetramer BETWEEN 2000 and 3000  UNION ALL
 SELECT tetramer, length(genomes), 75 as source_table FROM `PF05221.17_tetras` WHERE tetramer BETWEEN 2000 and 3000  UNION ALL
 SELECT tetramer, length(genomes), 76 as source_table FROM `PF06026.14_tetras` WHERE tetramer BETWEEN 2000 and 3000  UNION ALL
 SELECT tetramer, length(genomes), 77 as source_table FROM `PF13393.6_tetras` WHERE tetramer BETWEEN 2000 and 3000  UNION ALL
 SELECT tetramer, length(genomes), 78 as source_table FROM `PF17136.4_tetras` WHERE tetramer BETWEEN 2000 and 3000  UNION ALL
 SELECT tetramer, length(genomes), 79 as source_table FROM `PF01139.17_tetras` WHERE tetramer BETWEEN 2000 and 3000  ORDER BY tetramer, source_table

