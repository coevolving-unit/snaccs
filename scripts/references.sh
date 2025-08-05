# Allen M1

allen=https://idk-etl-prod-download-bucket.s3.amazonaws.com/aibs_human_m1_10x/matrix.csv
allen_meta=https://idk-etl-prod-download-bucket.s3.amazonaws.com/aibs_human_m1_10x/metadata.csv

wget $allen
wget $allen_meta

mv matrix.csv allen_M1_matrix.csv
mv metadata.csv allen_M1_metadata.csv

# Allen MCA

allen=https://idk-etl-prod-download-bucket.s3.amazonaws.com/aibs_human_ctx_smart-seq/matrix.csv
allen_meta=https://idk-etl-prod-download-bucket.s3.amazonaws.com/aibs_human_ctx_smart-seq/metadata.csv

wget $allen
wget $allen_meta

mv matrix.csv allen_MCA_matrix.csv
mv metadata.csv allen_MCA_metadata.csv

# BRAIN

# L23 IT (done)
ref=https://datasets.cellxgene.cziscience.com/2ef6cf58-f97b-489b-94f3-7ca6f6cdaf4d.rds
wget $ref
mv 2ef6cf58-f97b-489b-94f3-7ca6f6cdaf4d.rds BRAIN-L23IT.rds

# oligo (done)
ref=https://datasets.cellxgene.cziscience.com/5b25a593-1533-4a0a-84dc-7ca4db4d08a3.rds
wget $ref
mv 5b25a593-1533-4a0a-84dc-7ca4db4d08a3.rds BRAIN-oligo.rds

# L4IT (done)
ref=https://datasets.cellxgene.cziscience.com/5b7237a8-bbf5-480f-b0e0-8271a971b2e7.rds
wget $ref
mv 5b7237a8-bbf5-480f-b0e0-8271a971b2e7.rds BRAIN-L4IT.rds

# L5IT (done)
ref=https://datasets.cellxgene.cziscience.com/792816ee-7285-40c1-bb02-54b0a6f16e8a.rds
wget $ref
mv 792816ee-7285-40c1-bb02-54b0a6f16e8a.rds BRAIN-L5IT.rds

# SST (done)
ref=https://datasets.cellxgene.cziscience.com/6c81e7ca-009b-4157-a69b-3b6c6305ede3.rds
wget $ref
mv 6c81e7ca-009b-4157-a69b-3b6c6305ede3.rds BRAIN-SST.rds

# VIP (done)
ref=https://datasets.cellxgene.cziscience.com/d4df0234-1861-4440-bb2b-b6b7ea7e1197.rds
wget $ref
mv d4df0234-1861-4440-bb2b-b6b7ea7e1197.rds BRAIN-VIP.rds

# PVALB (done)
ref=https://datasets.cellxgene.cziscience.com/4f093a7a-6e4c-49b1-98ab-1d1d10acf988.rds
wget $ref
mv 4f093a7a-6e4c-49b1-98ab-1d1d10acf988.rds BRAIN-PVALB.rds

# L6IT (done)
ref=https://datasets.cellxgene.cziscience.com/125be0f9-84a9-4581-882d-aa65bbee6962.rds
wget $ref
mv 125be0f9-84a9-4581-882d-aa65bbee6962.rds BRAIN-L6IT.rds

# Astro (done)
ref=https://datasets.cellxgene.cziscience.com/5681afeb-e51d-4f87-b3df-611b9c2aabe8.rds
wget $ref
mv 5681afeb-e51d-4f87-b3df-611b9c2aabe8.rds BRAIN-Astro.rds

# OPC (done)
ref=https://datasets.cellxgene.cziscience.com/8ef6438e-9464-4bf8-85a8-2a577bd91506.rds
wget $ref
mv 8ef6438e-9464-4bf8-85a8-2a577bd91506.rds BRAIN-OPC.rds

# L6b (done)
ref=https://datasets.cellxgene.cziscience.com/5499e306-131a-4ada-8d74-ef53a86fe4a1.rds
wget $ref
mv 5499e306-131a-4ada-8d74-ef53a86fe4a1.rds BRAIN-L6b.rds

# LAMP5 (done)
ref=https://datasets.cellxgene.cziscience.com/caebc0eb-cea7-4b95-8123-be14111639f3.rds
wget $ref
mv caebc0eb-cea7-4b95-8123-be14111639f3.rds BRAIN-LAMP5.rds

# L6CT (done)
ref=https://datasets.cellxgene.cziscience.com/ffc75a9f-32f0-428e-bfa5-fd4e11017401.rds
wget $ref
mv ffc75a9f-32f0-428e-bfa5-fd4e11017401.rds BRAIN-L6CT.rds

# L56NP (done)
ref=https://datasets.cellxgene.cziscience.com/55262f6f-a2e8-4e92-b50d-840037ca979d.rds
wget $ref
mv 55262f6f-a2e8-4e92-b50d-840037ca979d.rds BRAIN-L56NP.rds

# L6ITCar3 (done)
ref=https://datasets.cellxgene.cziscience.com/78f91eb7-a8be-49f1-83fb-7c931059171c.rds
wget $ref
mv 78f91eb7-a8be-49f1-83fb-7c931059171c.rds BRAIN-L6ITCar3.rds

# Lamp5Lhx6 (done)
ref=https://datasets.cellxgene.cziscience.com/fafdb117-acfb-4a0b-bf80-3fc816ab47a2.rds
wget $ref
mv fafdb117-acfb-4a0b-bf80-3fc816ab47a2.rds BRAIN-Lamp5Lhx6.rds

# Sncg (done)
ref=https://datasets.cellxgene.cziscience.com/a844b19a-ed56-44fc-89f2-d1f69e8ab185.rds
wget $ref
mv a844b19a-ed56-44fc-89f2-d1f69e8ab185.rds BRAIN-Sncg.rds

# MicroPVM (done)
ref=https://datasets.cellxgene.cziscience.com/010f17d9-915d-45e6-9045-ca915694b4ca.rds
wget $ref
mv 010f17d9-915d-45e6-9045-ca915694b4ca.rds BRAIN-MicroPVM.rds

# Chandelier (done)
ref=https://datasets.cellxgene.cziscience.com/a86be20b-9531-439e-80a5-33adf94d7715.rds
wget $ref
mv a86be20b-9531-439e-80a5-33adf94d7715.rds BRAIN-Chandelier.rds

# Pax6 (done)
ref=https://datasets.cellxgene.cziscience.com/d90b8bb3-3fde-42be-b23b-dba98a9b4331.rds
wget $ref
mv d90b8bb3-3fde-42be-b23b-dba98a9b4331.rds BRAIN-Pax6.rds

# VLMC (done)
ref=https://datasets.cellxgene.cziscience.com/f707bb52-b8f1-4111-8970-78b07b172b55.rds
wget $ref
mv f707bb52-b8f1-4111-8970-78b07b172b55.rds BRAIN-VLMC.rds

# Endo (done)
ref=https://datasets.cellxgene.cziscience.com/fa2e9ed8-1629-463b-bf17-91c4738f692d.rds
wget $ref
mv fa2e9ed8-1629-463b-bf17-91c4738f692d.rds BRAIN-Endo.rds

# L5ET (done)
ref=https://datasets.cellxgene.cziscience.com/46166028-52fa-4d36-a774-e19d71d634f6.rds
wget $ref
mv 46166028-52fa-4d36-a774-e19d71d634f6.rds BRAIN-L5ET.rds

# SstChodl (done)
ref=https://datasets.cellxgene.cziscience.com/879fb8da-3581-4ed6-bd1c-ecd83094d7a2.rds
wget $ref
mv 879fb8da-3581-4ed6-bd1c-ecd83094d7a2.rds BRAIN-SstChodl.rds

# COPs (from Siletti channel)
ref=https://datasets.cellxgene.cziscience.com/98ec22fd-c595-4d77-b7b2-94f0cecd6024.rds
wget $ref
mv 98ec22fd-c595-4d77-b7b2-94f0cecd6024.rds BRAIN-COPs.rds
