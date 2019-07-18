import pandas as pd
import hashlib
import os
def file_as_bytes(file):
    with file:
        return file.read()

def extractMatrix(dirname):
    '''
    return a dataframe of the miRNA matrix, each row is the miRNA counts for a file_id
    '''
    count = 0

    miRNA_data = []
   
    for idname in os.listdir(dirname):
        # list all the ids
        if idname.find("-") != -1:
            idpath = dirname +"/" + idname

            # all the files in each id directory
            for filename in os.listdir(idpath):
                # check the miRNA file
                if filename.find("-") != -1:

                    filepath = idpath + "/" + filename
                    df = pd.read_csv(filepath,sep="\t")
                    # columns = ["miRNA_ID", "read_count"]
                    if count ==0:
                        # get the miRNA_IDs
                        miRNA_IDs = df.miRNA_ID.values.tolist()

                    id_miRNA_read_counts = [idname] + df.read_count.values.tolist()
                    miRNA_data.append(id_miRNA_read_counts)

                    count +=1
                    # print (df)
    columns = ["file_id"] + miRNA_IDs
    df = pd.DataFrame(miRNA_data, columns=columns)
    return df

def extractbc(df):
    #sample type
    df['label'] = df['cases.0.samples.0.sample_type']
    df.loc[df['cases.0.samples.0.sample_type'].str.contains("Solid"), 'label'] = 0
    df.loc[df['cases.0.samples.0.sample_type'].str.contains("Primary Tumor"), 'label'] = 1
    df.loc[df['cases.0.samples.0.sample_type'].str.contains("Recurrent Tumor"), 'label'] = 1
    df.loc[df['cases.0.samples.0.sample_type'].str.contains("Peripheral"), 'label'] = 1
    df.loc[df['cases.0.samples.0.sample_type'].str.contains("Metastatic"), 'label'] = 1
    df.loc[df['cases.0.samples.0.sample_type'].str.contains("Bone"), 'label'] = 1
    df.loc[df['cases.0.samples.0.sample_type'].str.contains("New"), 'label'] = 1
    df.loc[df['cases.0.samples.0.sample_type'].str.contains("Control"), 'label'] = 1
    df.loc[df['cases.0.samples.0.sample_type'].str.contains("Cell"), 'label'] = 1
   
    normal_count = df.loc[df.label == 1].shape[0]
    tumor_count = df.loc[df.label == 2].shape[0]
    columns = ['file_id','label']
    return df[columns]

def extractLabel_race(df):
    #race
    df['race'] = df['cases.0.demographic.race']
    df.loc[df['cases.0.demographic.race'].str.match("not reported"),'race'] = 1
    df.loc[df['cases.0.demographic.race'].str.match("asian"), 'race'] = 2
    df.loc[df['cases.0.demographic.race'].str.match("black or african american"), 'race'] = 3
    df.loc[df['cases.0.demographic.race'].str.match("white"), 'race'] = 4
    df.loc[df['cases.0.demographic.race'].str.match("native hawaiian or other pacific islander"), 'race'] = 5
    df.loc[df['cases.0.demographic.race'].str.match("american indian or alaska native"), 'race'] = 6
    df.loc[df['cases.0.demographic.race'].str.match("unclassified"), 'race'] = 7
    df.loc[df['cases.0.demographic.race'].str.match("other"), 'race'] = 8

    not_reported = df.loc[df.race == 1].shape[0]
    asian = df.loc[df.race == 2].shape[0]
    black = df.loc[df.race == 3].shape[0]
    white = df.loc[df.race == 4].shape[0]
    native_hawiian = df.loc[df.race == 5].shape[0]
    american_indian = df.loc[df.race == 6].shape[0]
    blank = df.loc[df.race == 7].shape[0]
    other = df.loc[df.race == 8].shape[0]
    columns = ['file_id','race']
    return df[columns]

def extractLabel_ethnicity(df):
    #ethnicity
    df['ethnicity'] = df['cases.0.demographic.ethnicity']
    df.loc[df['cases.0.demographic.ethnicity'].str.match("not reported"), 'ethnicity'] = 1
    df.loc[df['cases.0.demographic.ethnicity'].str.match("not hispanic or latino"), 'ethnicity'] = 2
    df.loc[df['cases.0.demographic.ethnicity'].str.match("hispanic or latino"), 'ethnicity'] = 3
    df.loc[df['cases.0.demographic.ethnicity'].str.match("unclassified"), 'ethnicity'] = 4
    df.loc[df['cases.0.demographic.ethnicity'].str.match("Unknown"), 'ethnicity'] = 5
   
    not_reported = df.loc[df.ethnicity == 1].shape[0]
    not_hispanic = df.loc[df.ethnicity == 2].shape[0]
    hispanic = df.loc[df.ethnicity == 3].shape[0]
    blank = df.loc[df.ethnicity == 4].shape[0]
    unknown = df.loc[df.ethnicity == 5].shape[0]
    columns = ['file_id','ethnicity']
    return df[columns]

def extractLabel_gender(df):
    #gender
    df['gender'] = df['cases.0.demographic.gender']
    df.loc[df['cases.0.demographic.gender'].str.contains("male"), 'gender'] = 1
    df.loc[df['cases.0.demographic.gender'].str.contains("female"), 'gender'] = 2
    df.loc[df['cases.0.demographic.gender'].str.contains("unclassified"), 'gender'] = 3
   
    male = df.loc[df.gender == 1].shape[0]
    female = df.loc[df.gender == 2].shape[0]
    blank = df.loc[df.gender == 3].shape[0]
    columns = ['file_id','gender']
    return df[columns]

def extractLabel_age(df):
    #patient age
    df['age'] = df['cases.0.diagnoses.0.age_at_diagnosis']
    df.loc[((df['cases.0.diagnoses.0.age_at_diagnosis']/365)>0) & ((df['cases.0.diagnoses.0.age_at_diagnosis']/365)<=10), 'age'] = 1
    df.loc[((df['cases.0.diagnoses.0.age_at_diagnosis']/365)>10) & ((df['cases.0.diagnoses.0.age_at_diagnosis']/365)<=20), 'age'] = 2
    df.loc[((df['cases.0.diagnoses.0.age_at_diagnosis']/365)>20) & ((df['cases.0.diagnoses.0.age_at_diagnosis']/365)<=30), 'age'] = 3
    df.loc[((df['cases.0.diagnoses.0.age_at_diagnosis']/365)>30) & ((df['cases.0.diagnoses.0.age_at_diagnosis']/365)<=40), 'age'] = 4
    df.loc[((df['cases.0.diagnoses.0.age_at_diagnosis']/365)>40) & ((df['cases.0.diagnoses.0.age_at_diagnosis']/365)<=50), 'age'] = 5
    df.loc[((df['cases.0.diagnoses.0.age_at_diagnosis']/365)>50) & ((df['cases.0.diagnoses.0.age_at_diagnosis']/365)<=60), 'age'] = 6
    df.loc[((df['cases.0.diagnoses.0.age_at_diagnosis']/365)>60) & ((df['cases.0.diagnoses.0.age_at_diagnosis']/365)<=70), 'age'] = 7
    df.loc[((df['cases.0.diagnoses.0.age_at_diagnosis']/365)>70) & ((df['cases.0.diagnoses.0.age_at_diagnosis']/365)<=80), 'age'] = 8
    df.loc[((df['cases.0.diagnoses.0.age_at_diagnosis']/365)>80) & ((df['cases.0.diagnoses.0.age_at_diagnosis']/365)<=90), 'age'] = 9
    df.loc[((df['cases.0.diagnoses.0.age_at_diagnosis']/365)>90) & ((df['cases.0.diagnoses.0.age_at_diagnosis']/365)<=100), 'age'] = 10
    df.loc[df['cases.0.diagnoses.0.age_at_diagnosis'] == 0, 'age'] = 11
   
    upto10 = df.loc[df.age == 1].shape[0]
    upto20 = df.loc[df.age == 2].shape[0]
    upto30 = df.loc[df.age == 3].shape[0]
    upto40 = df.loc[df.age == 4].shape[0]
    upto50 = df.loc[df.age == 5].shape[0]
    upto60 = df.loc[df.age == 6].shape[0]
    upto70 = df.loc[df.age == 7].shape[0]
    upto80 = df.loc[df.age == 8].shape[0]
    upto90 = df.loc[df.age == 9].shape[0]
    upto100 = df.loc[df.age == 10].shape[0]
    unclassified = df.loc[df.age == 11].shape[0]
    columns = ['file_id','age']
    return df[columns]

def extractSiteLabel(df):
        #cancer type
	df['cancer_type'] = df['cases.0.project.disease_type']   
	df.loc[df['cases.0.project.disease_type'].str.contains("Breast"), 'cancer_type'] = 0
	df.loc[df['cases.0.project.disease_type'].str.contains("Lung"), 'cancer_type'] = 1
	df.loc[df['cases.0.project.disease_type'].str.contains("Kidney"), 'cancer_type'] = 2
	df.loc[df['cases.0.project.disease_type'].str.contains("Ovarian"), 'cancer_type'] = 3
	df.loc[df['cases.0.project.disease_type'].str.contains("Brain"), 'cancer_type'] = 4
	df.loc[df['cases.0.project.disease_type'].str.contains("Stomach"), 'cancer_type'] = 5
	df.loc[df['cases.0.project.disease_type'].str.contains("Head and Neck"), 'cancer_type'] = 6
	df.loc[df['cases.0.project.disease_type'].str.contains("Bladder"), 'cancer_type'] = 7
	df.loc[df['cases.0.project.disease_type'].str.contains("Skin"), 'cancer_type'] = 8
	df.loc[df['cases.0.project.disease_type'].str.contains("Leukemia"), 'cancer_type'] = 9
	df.loc[df['cases.0.project.disease_type'].str.contains("Cervical"), 'cancer_type'] = 10
	df.loc[df['cases.0.project.disease_type'].str.contains("Liver"), 'cancer_type'] = 11
	df.loc[df['cases.0.project.disease_type'].str.contains("Uterine"), 'cancer_type'] = 12
	df.loc[df['cases.0.project.disease_type'].str.contains("Colon"), 'cancer_type'] = 13
	df.loc[df['cases.0.project.disease_type'].str.contains("Prostate"), 'cancer_type'] = 14
	df.loc[df['cases.0.project.disease_type'].str.contains("Thyroid"), 'cancer_type'] = 15
	df.loc[df['cases.0.project.disease_type'].str.contains("Rhabdoid"), 'cancer_type'] = 16
	df.loc[df['cases.0.project.disease_type'].str.contains("Thymoma"), 'cancer_type'] = 17
	df.loc[df['cases.0.project.disease_type'].str.contains("Cholangiocarcinoma"), 'cancer_type'] = 18
	df.loc[df['cases.0.project.disease_type'].str.contains("Uveal"), 'cancer_type'] = 19
	df.loc[df['cases.0.project.disease_type'].str.contains("Neoplasm"), 'cancer_type'] = 20
	df.loc[df['cases.0.project.disease_type'].str.contains("Glioblastoma"), 'cancer_type'] = 21
	df.loc[df['cases.0.project.disease_type'].str.contains("Sarcoma"), 'cancer_type'] = 22
	df.loc[df['cases.0.project.disease_type'].str.contains("Testicular"), 'cancer_type'] = 23
	df.loc[df['cases.0.project.disease_type'].str.contains("Rectum"), 'cancer_type'] = 24
	df.loc[df['cases.0.project.disease_type'].str.contains("Pancreatic"), 'cancer_type'] = 25
	df.loc[df['cases.0.project.disease_type'].str.contains("Esophageal"), 'cancer_type'] = 26
	df.loc[df['cases.0.project.disease_type'].str.contains("Pheochromocytoma"), 'cancer_type'] = 27
	df.loc[df['cases.0.project.disease_type'].str.contains("Adrenocortical"), 'cancer_type'] = 27
	df.loc[df['cases.0.project.disease_type'].str.contains("Mesothelioma"), 'cancer_type'] = 28
	df.loc[df['cases.0.project.disease_type'].str.contains("Wilms"), 'cancer_type'] = 29

	breast_count = df.loc[df.cancer_type == 0].shape[0]
	lung_count = df.loc[df.cancer_type == 1].shape[0]
	kidney_count = df.loc[df.cancer_type == 2].shape[0]
	ovary_count = df.loc[df.cancer_type == 3].shape[0]
	brain_count = df.loc[df.cancer_type == 4].shape[0]
	stom_count = df.loc[df.cancer_type == 5].shape[0]
	hn_count = df.loc[df.cancer_type == 6].shape[0]
	blad_count = df.loc[df.cancer_type == 7].shape[0]
	skin_count = df.loc[df.cancer_type == 8].shape[0]
	leuk_count = df.loc[df.cancer_type == 9].shape[0]
	sar_count = df.loc[df.cancer_type == 10].shape[0]
	cervix_count = df.loc[df.cancer_type == 11].shape[0]
	liver_count = df.loc[df.cancer_type == 12].shape[0]
	testis_count = df.loc[df.cancer_type == 13].shape[0]
	uterus_count = df.loc[df.cancer_type == 14].shape[0]
	colon_count = df.loc[df.cancer_type == 15].shape[0]
	rectum_count = df.loc[df.cancer_type == 16].shape[0]
	panc_count = df.loc[df.cancer_type == 17].shape[0]
	esop_count = df.loc[df.cancer_type == 18].shape[0]
	adg_count = df.loc[df.cancer_type == 19].shape[0]
	pluera_count = df.loc[df.cancer_type == 20].shape[0]
	prostate_count = df.loc[df.cancer_type == 21].shape[0]
	wilms_count = df.loc[df.cancer_type == 22].shape[0]
	thyroid_count = df.loc[df.cancer_type == 23].shape[0]
	rhab_count = df.loc[df.cancer_type == 24].shape[0]
	thy_count = df.loc[df.cancer_type == 25].shape[0]
	bile_count = df.loc[df.cancer_type == 26].shape[0]
	uveal_count = df.loc[df.cancer_type == 27].shape[0]
	lymph_count = df.loc[df.cancer_type == 28].shape[0]
	glio_count = df.loc[df.cancer_type == 29].shape[0]
	columns = ['file_id','cancer_type']
	return df[columns]

if __name__ == '__main__':

    data_dir ="/home/ubuntu/cancer/"
    # Input directory and label file. The directory that holds the data. Modify this when use.
    label_file = data_dir + "files_meta.tsv"
    outputfile = "/home/ubuntu/cancer/miRNA_matrix.csv"
    df = pd.read_csv(label_file, sep="\t")

    #fill NaN values in the dataframe
    df[['cases.0.demographic.ethnicity']] = df[['cases.0.demographic.ethnicity']].fillna("unclassified")
    df[['cases.0.demographic.gender']] = df[['cases.0.demographic.gender']].fillna("unclassified")
    df[['cases.0.demographic.race']] = df[['cases.0.demographic.race']].fillna("unclassified")
    df[['cases.0.diagnoses.0.age_at_diagnosis']] = df[['cases.0.diagnoses.0.age_at_diagnosis']].fillna(0)
    
    # extract data
    matrix_df = extractMatrix("/home/ubuntu/cancer/live_miRNA")
    label_df = extractbc(df)
    label1_df = extractLabel_race(df)
    label2_df = extractLabel_ethnicity(df)
    label3_df = extractLabel_gender(df)
    label4_df = extractLabel_age(df)
    label5_df = extractSiteLabel(df)
    lst=["Lung","Breast","Kidney","Stomach","Head and Neck","Prostate","Liver","Thyroid"]

    #computing csv files for different cancer types
    for i in range(len(lst)):
        print("Computing matrix for: ",lst[i])
        temp_matrix2=df[df['cases.0.project.disease_type'].str.contains(lst[i])]
    
        lc_matrix=temp_matrix2[['file_id']].copy()
    
        outputfile1 = "/home/ubuntu/cancer/"+lst[i]+"_matrix.csv"
        result1 = pd.merge(pd.merge(pd.merge(pd.merge(pd.merge(pd.merge(pd.merge(lc_matrix, matrix_df, on='file_id', how="left"), label5_df, on='file_id', how="left"), label4_df, on='file_id', how="left"), label3_df, on='file_id', how="left"), label2_df, on='file_id', how="left"), label1_df, on='file_id', how="left"), label_df, on='file_id', how="left")
        result1.to_csv(outputfile1, index=False)

    result = pd.merge(pd.merge(pd.merge(pd.merge(pd.merge(pd.merge(matrix_df, label5_df, on='file_id', how="left"), label4_df, on='file_id', how="left"), label3_df, on='file_id', how="left"), label2_df, on='file_id', how="left"), label1_df, on='file_id', how="left"), label_df, on='file_id', how="left")
    result.to_csv(outputfile, index=False)
