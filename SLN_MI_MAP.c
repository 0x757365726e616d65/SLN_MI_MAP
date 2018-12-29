//
//  SLN_MI_MAP.c / SLN_MI_MAP.h  //
//
//  Created by MingleZhao on 2018/12/21.
//  Copyright © 2018 Peking University. All rights reserved.
//

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <complex.h>
#include <sys/types.h>
#include <dirent.h>

#define M2L_16(x) ((unsigned short) \
((((unsigned short) ((x) & 0x00ff)) << 8) | \
(((unsigned short) ((x) & 0xff00)) >> 8)))

#define M2L_32(x) ((unsigned short) \
((((unsigned short) ((x) & 0x000000ff)) << 24) | \
(((unsigned short) ((x) & 0x0000ff00)) << 8) | \
(((unsigned short) ((x) & 0x00ff0000)) >> 8) | \
(((unsigned short) ((x) & 0xff000000)) >> 24)))
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int SLN_MI_MAP(char *DataPath, char *OutputPath);  // char *CalMode //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int SLN_MI_MAP(char *DataPath, char *OutputPath)  // char *CalMode //
{
    //__Directory_Path_End_Check__//
    unsigned short DataPathLen = strlen(DataPath);
    unsigned short OutputPathLen = strlen(OutputPath);
    
    char *NewDataPath = (char *) calloc(512, sizeof(char));
    char *NewDataPath0 = (char *) calloc(512, sizeof(char));
    char *NewOutputPath = (char *) calloc(512, sizeof(char));
    char *NewOutputPath0 = (char *) calloc(512, sizeof(char));
    if(NewDataPath == NULL || NewDataPath0 == NULL || NewOutputPath == NULL || NewOutputPath0 == NULL)
    {
        printf("\nError in allocating memory. \n");
        return -1;
    }
    
    strcpy(NewDataPath, DataPath);
    if(*((char *) (NewDataPath + DataPathLen - 1)) != '/')
    {
        *((char *) (NewDataPath + DataPathLen)) = '/';
        *((char *) (NewDataPath + DataPathLen + 1)) = '\0';
    }
    strcpy(NewOutputPath, OutputPath);
    if(*((char *) (NewOutputPath + OutputPathLen - 1)) != '/')
    {
        *((char *) (NewOutputPath + OutputPathLen)) = '/';
        *((char *) (NewOutputPath + OutputPathLen + 1)) = '\0';
    }
    strcpy(NewDataPath0, NewDataPath);
    strcpy(NewOutputPath0, NewOutputPath);
    //
    
    //__Definations_of_Pointers__//
    FILE *fd_Data = NULL, *fd_Label = NULL, *fd_Out = NULL;
    struct dirent *DataDirInfo = NULL;
    DIR *DataDP = opendir(NewDataPath);
    if(DataDP == NULL)
    {
        printf("\nError in opening directory. => 1ST \n");
        return -1;
    }
    char *DataName = (char *) calloc(512, sizeof(char));
    char *LabelName = (char *) calloc(512, sizeof(char));
    char **DataNameList = (char **) calloc(1 * sizeof(char *));
    if((DataName == NULL) || (LabelName == NULL) || (DataNameList == NULL))
    {
        printf("\nError in allocating memory. \n");
        return -1;
    }
    //
    
    short Flag = 0;
    unsigned int DataNum = 0, Nth = 0;
    //__Counting_the_Number_of_Data__//
    //__Writing_the_Data_Name_List__//
    while((DataDirInfo = readdir(DataDP)))
    {
        if(((strstr(DataDirInfo -> d_name, ".img")) || (strstr(DataDirInfo -> d_name, ".IMG"))) && (strstr(DataDirInfo -> d_name, "MI_MAP")))
        {
            DataNum += 1;
            strcpy(DataName, DataDirInfo -> d_name);
            
            if(!(DataNameList = (char **) realloc(DataNameList, DataNum * sizeof(char *))))
            {
                printf("\nError in reallocating memory. \n");
                return -1;
            }
            
            if(!(*(DataNameList + DataNum - 1) = (char *) calloc(strlen(DataName) + 1, sizeof(char))))
            {
                printf("\nError in allocating memory. \n");
                return -1;
            }
            
            strcpy(*(DataNameList + DataNum - 1), DataName);
        }
    }
    free(DataDirInfo);
    closedir(DataDP);
    DataDirInfo = NULL;
    DataDP = NULL;
    //
    //
    
    //__Reading_Data_and_Label_Files__//
    for(Nth = 1; Nth <= DataNum; Nth += 1)
    {
        printf("\nSum = %d, Now = %d\n", DataNum, Nth);
        
        //__Declarations_of_Variables_in_Reading_Data_and_Label_Files__//
        long ptr_shift = 0, start = 0, start_geo = 0;
        char *p_inline = NULL;
        //__Declarations_of_Label_Parameters__//
        int lines = 0, samples = 0, bytes = 0;
        double latnst = 0, lonwst = 0;
        double map_resolution = 0, scaling_factor = 0, offset = 0;
        double line_projection_offset = 0, sample_projection_offset = 0;
        char line[1024] = {'\0'};
        //
        //
        
        strcpy(DataName, *(DataNameList + Nth - 1));
        memset(NewDataPath0, 0, 512);
        strcpy(NewDataPath0, NewDataPath);
        if((fd_Data = fopen(strcat(NewDataPath0, DataName), "r")) == NULL)
        {
            printf("\nError in opening data. => Data = %s\n", DataName);
            return -1;
        }
        printf("Data:\t%s\n", DataName);
        
        //__Judging_Data_in_New_or_Previous_Format//
        fseek(fd_Data, 0, SEEK_SET);
        if(strstr(fgets(line, sizeof(line), fd_Data), "PDS_VERSION_ID"))
        {
            Flag = 0;
            fd_Label = fd_Data;
            strcpy(LabelName, DataName);
        }
        else
        {
            Flag = 1;
            strcpy(LabelName, DataName);
            LabelName[strlen(LabelName) - 3] = 'l';
            LabelName[strlen(LabelName) - 2] = 'b';
            LabelName[strlen(LabelName) - 1] = 'l';
            
            memset(NewDataPath0, 0, 512);
            strcpy(NewDataPath0, NewDataPath);
            if((fd_Label = fopen(strcat(NewDataPath0, LabelName), "r")) == NULL)
            {
                printf("\nError in opening the label file. => Label = %s\n", LabelName);
                return -1;
            }
        }
        printf("Label:\t%s\n", LabelName);
    
        //__Label_Parameters__//
        short Flag_SampleType = 0;
        fseek(fd_Label, 0, SEEK_SET);
        while((strcmp(fgets(line, sizeof(line), fd_Label), "END\r\n")))
        {
            if((strstr(line, "MSB_INTEGER")))
            {
                Flag_SampleType = 1;
            }
            
            if((strstr(line, "SAMPLE_BITS")) && (Flag_SampleType == 1))
            {
                p_inline = strchr(line, '=');
                sscanf(p_inline + 1, "%d", &bytes);
                printf("SAMPLE_BITS\t=\t%d\n", bytes);
                bytes /= 8;
            }
            
            if((strstr(line, "SCALING_FACTOR")))
            {
                p_inline = strchr(line, '=');
                sscanf(p_inline + 1, "%lf", &scaling_factor);
                printf("SCALING_FACTOR\t=\t%lf\n", scaling_factor);
            }
            
            if((strstr(line, "OFFSET")) && (!strstr(line, "_")) && (!strstr(line, "*")))
            {
                p_inline = strchr(line, '=');
                sscanf(p_inline + 1, "%lf", &offset);
                printf("OFFSET\t=\t%lf\n", offset);
            }
            
            if((strstr(line, "LINE_LAST_PIXEL")))
            {
                p_inline = strchr(line, '=');
                sscanf(p_inline + 1, "%d", &lines);
                printf("LINE_LAST_PIXEL\t=\t%d\n", lines);
            }
            
            if((strstr(line, "SAMPLE_LAST_PIXEL")))
            {
                p_inline = strchr(line, '=');
                sscanf(p_inline + 1, "%d", &samples);
                printf("SAMPLE_LAST_PIXEL\t=\t%d\n", samples);
            }
            
            if((strstr(line, "MAP_RESOLUTION")))
            {
                p_inline = strchr(line, '=');
                sscanf(p_inline + 1, "%lf", &map_resolution);
                printf("MAP_RESOLUTION\t=\t%lf\n", map_resolution);
            }
            
            if((strstr(line, "MAXIMUM_LATITUDE")))
            {
                p_inline = strchr(line, '=');
                sscanf(p_inline + 1, "%lf", &latnst);
                printf("MAXIMUM_LATITUDE\t=\t%lf\n", latnst);
            }
            
            if((strstr(line, "WESTERNMOST_LONGITUDE")))
            {
                p_inline = strchr(line, '=');
                sscanf(p_inline + 1, "%lf", &lonwst);
                printf("WESTERNMOST_LONGITUDE\t=\t%lf\n", lonwst);
            }
            
            if((strstr(line, "LINE_PROJECTION_OFFSET")))
            {
                p_inline = strchr(line, '=');
                sscanf(p_inline + 1, "%lf", &line_projection_offset);
                printf("LINE_PROJECTION_OFFSET\t=\t%lf\n", line_projection_offset);
            }
            
            if((strstr(line, "SAMPLE_PROJECTION_OFFSET")))
            {
                p_inline = strchr(line, '=');
                sscanf(p_inline + 1, "%lf", &sample_projection_offset);
                printf("SAMPLE_PROJECTION_OFFSET\t=\t%lf\n", sample_projection_offset);
            }
            
            if((strstr(line, "^GEOMETRIC_DATA_ALTITUDE")))
            {
                if(Flag == 0)
                    p_inline = strchr(line, '=');
                else if(Flag == 1)
                    p_inline = strchr(line, ',');
                sscanf(p_inline + 1, "%ld", &start_geo);
                printf("^GEOMETRIC_DATA_ALTITUDE\t=\t%ld\n", start_geo);
            }
            
            if((strstr(line, "^IMAGE")))
            {
                if(Flag == 0)
                    p_inline = strchr(line, '=');
                else if(Flag == 1)
                    p_inline = strchr(line, ',');
                sscanf(p_inline + 1, "%ld", &start);
                printf("^IMAGE\t=\t%ld\n", start);
            }
            
            memset(line, 0, sizeof(line));
        }
        //
        
        //__Decalarations_of_Parameters__//
        //__FeOTiO2_Algorithm__//
        double x0fe = 0.04, y0fe = 1.39;  // Lemelin et al. 2015 //
        double x0ti = -0.108, y0ti = 0.208;  // Otake et al. 2012 //
        double x0 = 0.037, y0 = 1.250;  //  //
        //
        
        //__Calculation_Parameters__//
        int i = 0, j = 0, BandNum = 0;
        short value = 0;
        double lat = 0, lon = 0;
        double theta_fe1 = 0, theta_fe2 = 0, theta_ti = 0;
        double lgtan = 0, b0 = 0;
        double _Complex e = 0, e0 = 0, A = 0;
        //
        
        //__Multiband_Reflectance_Matrix__//
        double **r415 = (double **) calloc(lines, sizeof(double *));
        if(r415 == NULL)
        {
            printf("Error in allocating memory. \n");
            return -1;
        }
        for(i = 0; i < lines; i += 1)
        {
            r415[i] = (double *) calloc(samples, sizeof(double));
            if(r415[i] == NULL)
            {
                printf("Error in allocating memory. \n");
                return -1;
            }
        }
        
        double **r750 = (double **) calloc(lines, sizeof(double *));
        if(r750 == NULL)
        {
            printf("Error in allocating memory. \n");
            return -1;
        }
        for(i = 0; i < lines; i += 1)
        {
            r750[i] = (double *) calloc(samples, sizeof(double));
            if(r750[i] == NULL)
            {
                printf("Error in allocating memory. \n");
                return -1;
            }
        }
        
        double **r900 = (double **) calloc(lines, sizeof(double *));
        if(r900 == NULL)
        {
            printf("Error in allocating memory. \n");
            return -1;
        }
        for(i = 0; i < lines; i += 1)
        {
            r900[i] = (double *) calloc(samples, sizeof(double));
            if(r900[i] == NULL)
            {
                printf("Error in allocating memory. \n");
                return -1;
            }
        }
        
        double **r950 = (double **) calloc(lines, sizeof(double *));
        if(r950 == NULL)
        {
            printf("Error in allocating memory. \n");
            return -1;
        }
        for(i = 0; i < lines; i += 1)
        {
            r950[i] = (double *) calloc(samples, sizeof(double));
            if(r950[i] == NULL)
            {
                printf("Error in allocating memory. \n");
                return -1;
            }
        }
        
        double **r1000 = (double **) calloc(lines, sizeof(double *));
        if(r1000 == NULL)
        {
            printf("Error in allocating memory. \n");
            return -1;
        }
        for(i = 0; i < lines; i += 1)
        {
            r1000[i] = (double *) calloc(samples, sizeof(double));
            if(r1000[i] == NULL)
            {
                printf("Error in allocating memory. \n");
                return -1;
            }
        }
        //
        
        //__Output_Products_Matrix__//
        double **fe = (double **) calloc(lines, sizeof(double *));
        if(fe == NULL)
        {
            printf("Error in allocating memory. \n");
            return -1;
        }
        for(i = 0; i < lines; i += 1)
        {
            fe[i] = (double *) calloc(samples, sizeof(double));
            if(fe[i] == NULL)
            {
                printf("Error in allocating memory. \n");
                return -1;
            }
        }
        
        double **ti = (double **) calloc(lines, sizeof(double *));
        if(ti == NULL)
        {
            printf("Error in allocating memory. \n");
            return -1;
        }
        for(i = 0; i < lines; i += 1)
        {
            ti[i] = (double *) calloc(samples, sizeof(double));
            if(ti[i] == NULL)
            {
                printf("Error in allocating memory. \n");
                return -1;
            }
        }
        
        double **den = (double **) calloc(lines, sizeof(double *));
        if(den == NULL)
        {
            printf("Error in allocating memory. \n");
            return -1;
        }
        for(i = 0; i < lines; i += 1)
        {
            den[i] = (double *) calloc(samples, sizeof(double));
            if(den[i] == NULL)
            {
                printf("Error in allocating memory. \n");
                return -1;
            }
        }
        
        double **re_diel = (double **) calloc(lines, sizeof(double *));
        if(re_diel == NULL)
        {
            printf("Error in allocating memory. \n");
            return -1;
        }
        for(i = 0; i < lines; i += 1)
        {
            re_diel[i] = (double *) calloc(samples, sizeof(double));
            if(re_diel[i] == NULL)
            {
                printf("Error in allocating memory. \n");
                return -1;
            }
        }
        
        double **im_diel = (double **) calloc(lines, sizeof(double *));
        if(im_diel == NULL)
        {
            printf("Error in allocating memory. \n");
            return -1;
        }
        for(i = 0; i < lines; i += 1)
        {
            im_diel[i] = (double *) calloc(samples, sizeof(double));
            if(im_diel[i] == NULL)
            {
                printf("Error in allocating memory. \n");
                return -1;
            }
        }
        
        double **d0 = (double **) calloc(lines, sizeof(double *));
        if(d0 == NULL)
        {
            printf("Error in allocating memory. \n");
            return -1;
        }
        for(i = 0; i < lines; i += 1)
        {
            d0[i] = (double *) calloc(samples, sizeof(double));
            if(d0[i] == NULL)
            {
                printf("Error in allocating memory. \n");
                return -1;
            }
        }
        
        double **omat = (double **) calloc(lines, sizeof(double *));
        if(omat == NULL)
        {
            printf("Error in allocating memory. \n");
            return -1;
        }
        for(i = 0; i < lines; i += 1)
        {
            omat[i] = (double *) calloc(samples, sizeof(double));
            if(omat[i] == NULL)
            {
                printf("Error in allocating memory. \n");
                return -1;
            }
        }
        //
        //
        
        //__Read_Multiband_Reflectance__//
        fseek(fd_Data, 0, SEEK_SET);
        //__r415__//
        BandNum += 1;
        ptr_shift = (start - 1) + (lines * samples * bytes * (BandNum - 1));
        fseek(fd_Data, ptr_shift, SEEK_SET);
        
        for(i = 0; i < lines; i += 1)
        {
            for (j = 0; j < samples; j += 1)
            {
                fread(&value, bytes, 1, fd_Data);
                value = M2L_16(value);
                r415[i][j] = value * scaling_factor + offset;
            }
        }
        //
        
        //__r750__//
        BandNum += 1;
        ptr_shift = (start - 1) + (lines * samples * bytes * (BandNum - 1));
        fseek(fd_Data, ptr_shift, SEEK_SET);
        
        for(i = 0; i < lines; i += 1)
        {
            for (j = 0; j < samples; j += 1)
            {
                fread(&value, bytes, 1, fd_Data);
                value = M2L_16(value);
                r750[i][j] = value * scaling_factor + offset;
            }
        }
        //
        
        //__r900__//
        BandNum += 1;
        ptr_shift = (start - 1) + (lines * samples * bytes * (BandNum - 1));
        fseek(fd_Data, ptr_shift, SEEK_SET);
        
        for(i = 0; i < lines; i += 1)
        {
            for (j = 0; j < samples; j += 1)
            {
                fread(&value, bytes, 1, fd_Data);
                value = M2L_16(value);
                r900[i][j] = value * scaling_factor + offset;
            }
        }
        //
        
        //__r950__//
        BandNum += 1;
        ptr_shift = (start - 1) + (lines * samples * bytes * (BandNum - 1));
        fseek(fd_Data, ptr_shift, SEEK_SET);
        
        for(i = 0; i < lines; i += 1)
        {
            for (j = 0; j < samples; j += 1)
            {
                fread(&value, bytes, 1, fd_Data);
                value = M2L_16(value);
                r950[i][j] = value * scaling_factor + offset;
            }
        }
        //
        
        //__r1000__//
        BandNum += 1;
        ptr_shift = (start - 1) + (lines * samples * bytes * (BandNum - 1));
        fseek(fd_Data, ptr_shift, SEEK_SET);
        
        for(i = 0; i < lines; i += 1)
        {
            for (j = 0; j < samples; j += 1)
            {
                fread(&value, bytes, 1, fd_Data);
                value = M2L_16(value);
                r1000[i][j] = value * scaling_factor + offset;
            }
        }
        //
        //
        
        //__Computations_of_Output_Products__//
        //__FeO__//
        memset(NewOutputPath0, 0, 512);
        strcpy(NewOutputPath0, NewOutputPath);
        if((fd_Out = fopen(strcat(NewOutputPath0, "FeO.txt"), "a")) == NULL)
        {
            printf("\nError in opening files. \n");
            return -1;
        }
        
        for (i = 0; i < lines; i += 1)
        {
            lat = latnst -  ((double) i) / map_resolution - 0.5 / map_resolution;
            for (j = 0; j < samples; j += 1)
            {
                lon = lonwst + ((double) j) / map_resolution + 0.5 / map_resolution;
                theta_fe1 = (-1) * atan((r950[i][j] / r750[i][j] - y0fe) / (r750[i][j] - x0fe));
                theta_fe2 = 0.0656 * exp(3.6681 * theta_fe1);
                fe[i][j] = 1.0708 * theta_fe2 - 0.3986;
                
                if(fe[i][j] >= 0 && fe[i][j] <= 100)
                    fprintf(fd_Out, "%.8lf\t%.8lf\t%.8lf\n", lon, lat, fe[i][j]);
            }
        }
        fclose(fd_Out);
        //
        
        //__TiO2__//
        memset(NewOutputPath0, 0, 512);
        strcpy(NewOutputPath0, NewOutputPath);
        if((fd_Out = fopen(strcat(NewOutputPath0, "TiO2.txt"), "a")) == NULL)
        {
            printf("\nError in opening files. \n");
            return -1;
        }
        
        for (i = 0; i < lines; i += 1)
        {
            lat = latnst -  ((double) i) / map_resolution - 0.5 / map_resolution;
            for (j = 0; j < samples; j += 1)
            {
                lon = lonwst + ((double) j) / map_resolution + 0.5 / map_resolution;
                theta_ti = atan((r415[i][j] / r750[i][j] - y0ti) / (r750[i][j] - x0ti));
                ti[i][j] = 0.72 * pow(theta_ti, 14.964);
                
                if(ti[i][j] >= 0 && ti[i][j] <= 100)
                    fprintf(fd_Out, "%.8lf\t%.8lf\t%.8lf\n", lon, lat, ti[i][j]);
            }
        }
        fclose(fd_Out);
        //
        
        //__Density__//
        memset(NewOutputPath0, 0, 512);
        strcpy(NewOutputPath0, NewOutputPath);
        if((fd_Out = fopen(strcat(NewOutputPath0, "Density.txt"), "a")) == NULL)
        {
            printf("\nError in opening files. \n");
            return -1;
        }
        
        for (i = 0; i < lines; i += 1)
        {
            lat = latnst -  ((double) i) / map_resolution - 0.5 / map_resolution;
            for (j = 0; j < samples; j += 1)
            {
                lon = lonwst + ((double) j) / map_resolution + 0.5 / map_resolution;
                den[i][j] = 0.0273 * fe[i][j] + 0.011 * ti[i][j] + 2.773;
                
                fprintf(fd_Out, "%.8lf\t%.8lf\t%.8lf\n", lon, lat, den[i][j]);
            }
        }
        fclose(fd_Out);
        //
        
        //__Real_Dielectric__//
        memset(NewOutputPath0, 0, 512);
        strcpy(NewOutputPath0, NewOutputPath);
        if((fd_Out = fopen(strcat(NewOutputPath0, "RealDielectric.txt"), "a")) == NULL)
        {
            printf("\nError in opening files. \n");
            return -1;
        }
        
        for (i = 0; i < lines; i += 1)
        {
            lat = latnst -  ((double) i) / map_resolution - 0.5 / map_resolution;
            for (j = 0; j < samples; j += 1)
            {
                lon = lonwst + ((double) j) / map_resolution + 0.5 / map_resolution;
                lgtan = -2.395 + 0.064 * ti[i][j];
                b0 = 2.75 * (pow(10, lgtan));
                e0 = 2.75 + b0*I;
                A = (den[i][j] * 0.55 / 1.7) * ((e0 - 1) / (e0 + 2));  // porosity == 0.45 //
                e = (2 * A + 1) / (1 - A);
                re_diel[i][j] = creal(e);
                
                fprintf(fd_Out, "%.8lf\t%.8lf\t%.8lf\n", lon, lat, re_diel[i][j]);
            }
        }
        fclose(fd_Out);
        //
        
        //__Imaginary_Dielectric__//
        memset(NewOutputPath0, 0, 512);
        strcpy(NewOutputPath0, NewOutputPath);
        if((fd_Out = fopen(strcat(NewOutputPath0, "ImaginaryDielectric.txt"), "a")) == NULL)
        {
            printf("\nError in opening files. \n");
            return -1;
        }
        
        for (i = 0; i < lines; i += 1)
        {
            lat = latnst -  ((double) i) / map_resolution - 0.5 / map_resolution;
            for (j = 0; j < samples; j += 1)
            {
                lon = lonwst + ((double) j) / map_resolution + 0.5 / map_resolution;
                lgtan = -2.395 + 0.064 * ti[i][j];
                b0 = 2.75 * (pow(10, lgtan));
                e0 = 2.75 + b0*I;
                A = (den[i][j] * 0.55 / 1.7) * ((e0 - 1) / (e0 + 2));  // porosity == 0.45 //
                e = (2 * A + 1) / (1 - A);
                im_diel[i][j] = cimag(e);
                
                fprintf(fd_Out, "%.8lf\t%.8lf\t%.8lf\n", lon, lat, im_diel[i][j]);
            }
        }
        fclose(fd_Out);
        //
        
        //__d0_Normalized_Penetration_Depth__//
        memset(NewOutputPath0, 0, 512);
        strcpy(NewOutputPath0, NewOutputPath);
        if((fd_Out = fopen(strcat(NewOutputPath0, "d0_NormalizedPenetrationDepth.txt"), "a")) == NULL)
        {
            printf("\nError in opening files. \n");
            return -1;
        }
        
        for (i = 0; i < lines; i += 1)
        {
            lat = latnst -  ((double) i) / map_resolution - 0.5 / map_resolution;
            for (j = 0; j < samples; j += 1)
            {
                lon = lonwst + ((double) j) / map_resolution + 0.5 / map_resolution;
                d0[i][j] = (sqrt(re_diel[i][j])) / (2 * 3.141592653589793 * im_diel[i][j]);
                
                fprintf(fd_Out, "%.8lf\t%.8lf\t%.8lf\n", lon, lat, d0[i][j]);
            }
        }
        fclose(fd_Out);
        //
        
        //__OMAT__//
        memset(NewOutputPath0, 0, 512);
        strcpy(NewOutputPath0, NewOutputPath);
        if((fd_Out = fopen(strcat(NewOutputPath0, "OMAT.txt"), "a")) == NULL)
        {
            printf("\nError in opening files. \n");
            return -1;
        }
        
        for (i = 0; i < lines; i += 1)
        {
            lat = latnst -  ((double) i) / map_resolution - 0.5 / map_resolution;
            for (j = 0; j < samples; j += 1)
            {
                lon = lonwst + ((double) j) / map_resolution + 0.5 / map_resolution;
                omat[i][j] = sqrt(pow((r750[i][j] - x0), 2) + pow((r950[i][j] / r750[i][j] - y0), 2));
                
                fprintf(fd_Out, "%.8lf\t%.8lf\t%.8lf\n", lon, lat, omat[i][j]);
            }
        }
        fclose(fd_Out);
        //
        //
        
        //__Release_Pointers_and_Memory__//
        for(i = 0; i < lines; i += 1)
        {
            free(r415[i]);
        }
        free(r415);
        r415 = NULL;
        
        for(i = 0; i < lines; i += 1)
        {
            free(r750[i]);
        }
        free(r750);
        r750 = NULL;
        
        for(i = 0; i < lines; i += 1)
        {
            free(r900[i]);
        }
        free(r900);
        r900 = NULL;
        
        for(i = 0; i < lines; i += 1)
        {
            free(r950[i]);
        }
        free(r950);
        r950 = NULL;
        
        for(i = 0; i < lines; i += 1)
        {
            free(r1000[i]);
        }
        free(r1000);
        r1000 = NULL;
        
        for(i = 0; i < lines; i += 1)
        {
            free(fe[i]);
        }
        free(fe);
        fe = NULL;
        
        for(i = 0; i < lines; i += 1)
        {
            free(ti[i]);
        }
        free(ti);
        ti = NULL;
        
        for(i = 0; i < lines; i += 1)
        {
            free(den[i]);
        }
        free(den);
        den = NULL;
        
        for(i = 0; i < lines; i += 1)
        {
            free(re_diel[i]);
        }
        free(re_diel);
        re_diel = NULL;
        
        for(i = 0; i < lines; i += 1)
        {
            free(im_diel[i]);
        }
        free(im_diel);
        im_diel = NULL;
        
        for(i = 0; i < lines; i += 1)
        {
            free(d0[i]);
        }
        free(d0);
        d0 = NULL;
        
        for(i = 0; i < lines; i += 1)
        {
            free(omat[i]);
        }
        free(omat);
        omat = NULL;
        //
        
        printf("Finished = %d\n", Nth);
    }
    //
    
    //__Release_Pointers__//
    fclose(fd_Data);
    fclose(fd_Label);
    fclose(fd_Out);
    free(DataName);
    free(LabelName);
    
    fd_Data = NULL;
    fd_Label = NULL;
    fd_Out = NULL;
    DataName = NULL;
    LabelName = NULL;
    //
    
    //__Release_Pointers_in_Data_Name_List__//
    for(Nth = 1; Nth <= DataNum; Nth += 1)
    {
        free(*(DataNameList + Nth - 1));
    }
    free(DataNameList);
    DataNameList = NULL;
    //
    
    //__Release_Pointers__//
    free(NewDataPath);
    free(NewDataPath0);
    free(NewOutputPath);
    free(NewOutputPath0);
    
    NewDataPath = NULL;
    NewDataPath0 = NULL;
    NewOutputPath = NULL;
    NewOutputPath0 = NULL;
    //
    
    return((void) (printf("\nSLN_MI_MAP => Done. \n")), 0);
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
