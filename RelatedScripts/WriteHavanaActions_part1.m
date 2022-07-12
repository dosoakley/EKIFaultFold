function [] = WriteHavanaActions_part1(filename,strike_range,dip_range,anisotropy_angle,...
    asymmetry,range,ratio,power,length_height_ratio)
%WriteHavanaActions Write the HavanaActions xml file.
%   This is a script to write out a Havana actions xml file for importing
%   an RMS model, with specified values of the kinematic model parameters.

%The fault and horizon names are currently hard-coded. So are the folder 
%names. This is just for modifying kinematic parameters within this 
%specific geologic model.

file = fopen(filename,'w');
fprintf(file,'<havana>\n');
fprintf(file,'	<project-settings>\n');
fprintf(file,'		<io-settings>\n');
fprintf(file,'			<input-directory>.</input-directory>\n');
fprintf(file,'			<output-directory>.</output-directory>\n');
fprintf(file,'		</io-settings>\n');
fprintf(file,'		<log-settings>\n');
fprintf(file,'			<log-file>\n');
fprintf(file,'				<file>logfile1.dat</file>\n');
fprintf(file,'				<level>5</level>\n');
fprintf(file,'			</log-file>\n');
fprintf(file,'			<screen-level>1</screen-level>\n');
fprintf(file,'		</log-settings>\n');
fprintf(file,'	</project-settings>\n');
fprintf(file,'	<import-rms-fault-data>\n');
fprintf(file,'		<input-structural-model-directory>DeformedModel</input-structural-model-directory>\n');
fprintf(file,'		<output-havana-structural-model-directory>HavanaStructuralModel1</output-havana-structural-model-directory>\n');
fprintf(file,'		<debug-output-fault-info-directory> Debug </debug-output-fault-info-directory>\n');
fprintf(file,'		<displacement-variogram>\n');
fprintf(file,'			<type> Spherical </type>\n');
fprintf(file,'			<strike-range> %f </strike-range>\n',strike_range);
fprintf(file,'			<dip-range> %f </dip-range>\n',dip_range);
fprintf(file,'			<anisotropy-angle> %f </anisotropy-angle>\n',anisotropy_angle);
fprintf(file,'		</displacement-variogram>\n');
fprintf(file,'		<default-displacement-settings>\n');
fprintf(file,'			<asymmetry> %f </asymmetry>\n',asymmetry);
fprintf(file,'			<range> %f </range>\n',range);
fprintf(file,'			<max-displacement-length-relation>\n');
fprintf(file,'				<ratio> %f </ratio>\n',ratio);
fprintf(file,'				<power> %f </power>\n',power);
fprintf(file,'			</max-displacement-length-relation>\n');
fprintf(file,'			<length-height-ratio> %f </length-height-ratio>\n',length_height_ratio);
fprintf(file,'		</default-displacement-settings>\n');
fprintf(file,'	</import-rms-fault-data>\n');
fprintf(file,'</havana>');
fclose(file);

end

