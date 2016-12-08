%-----------------------------------------------------------------------------------------------------------------------
% octave m file
%
% Author:		Jan Stegenga
% Company:		INCAS3
% Date:			3/10/2014
% Copyright:	INCAS3 2016
% License:		GNU GPLv3.0
% Requirements:	32bit octave (i used windows 7).
%
% linear or quadratic stabilization of gamma spectra in the range of 0 to 3000 keV, using several approaches.
%	Full Spectrum Analysis:	
%		- brute force
%		- sqp
%	Peak finding: 
%		- with expected peak locations
%		- without expected peak locations
% 
% The goal is to generate spectra that have an identical energy axis and ths can be compared to eachother. 
% Generally, the energy axis of a multi channel analyser depends on temperature (linear scaling) and has an offset and 
% quadratic component due to instrumentation limitations. These can be corrected for by rescaling the energy-axis such 
% that (1) features with known energy are locted correctly (peak finding algorithms) or (2) the entire spectrum matches
% to a mixture of (monte carlo) simulated spectra (FSA). This version ignores live time / dead time and has many fixed
% factors throughout the code that should have been neatly tucked into a configuration file along with what they do.
%
% input:
%	csv file with accumulated multichannel spectra of 1 or 2 hours (i.e. XXX counts ), usually per column:
%		data - 1024 or 2048 integers
%		t	 - date AND time up to seconds
% 	
% output:
%	csv file with stabilized data in same organization
%	several plots
%-----------------------------------------------------------------------------------------------------------------------

clc;
close all;
clear all;
pkg load signal
pkg load optim
warning off

%-----------------------------------------------------------------------------------------------------------------------
%rebin a stretched spectrum so that the binwidth is equal over the whole energy range
%-----------------------------------------------------------------------------------------------------------------------
function newspectrum = rebin( energyaxis, spectrum )
	%--------------
	%This function outputs a spectrum from 0 to 3000 MeV in equal bins
	%given a spectrum and the energy-axis with unequal bin widths
	%example:
	%	stab_s  = rebin( polyval( [quadratic, linear, offset], 1:1024 ), spectrum );
	%--------------

	%build a cumulative spectrum
	cumdata 	   = cumsum( spectrum .* double( energyaxis >= 0 ) ); 
	%interpolate the cumulative spectrum at points given by energyaxis
	newspectrum    = interp1( energyaxis, cumdata, linspace( 0, 3000, 1*length(energyaxis)+1 ), 'linear', 'extrap' ); 
	%take the difference to revert to spectrum
	newspectrum    = diff( newspectrum( 1:1:end ) );
end
%-----------------------------------------------------------------------------------------------------------------------
%perform FSA using a number of reference spectra and sqp (Broyden–Fletcher–Goldfarb–Shanno) / fminsearch (Nelder-Mead) optimizers
%-----------------------------------------------------------------------------------------------------------------------
function [stab_s, p, act_con, error] = FSA_optimizers( spectrum, references, order, verbose )
	N 	= length( spectrum );
	LO 	= 100;	%assumes 1024 channels!!
	HI 	= 900;	%assumes 1024 channels!!
	s 	= spectrum';
	
	%give relative emphasis to errors using weighting and windowing vector. 
	%logarithmic weighting effectively means we're working in linear space from here on
	w 	= logspace( -3, 0, N );			
	wc 	= reshape( w.*double( ( 1:N > LO ) & ( 1:N < HI ) ), size( s' ) );  
	
	%apply w to references
    r = bsxfun( @times, references, wc );
	%energy axis
	e = linspace( 0, 3000, N );
	%inititial values for activity concentration in logspace:
	nref        = min( size( r ) );
	ncountsref  = sum( r(:) ) / nref;
	ncounts		= sum( s*wc );
	multiplier  = ncounts / ncountsref / nref;
	if order == 2
		p = [multiplier*ones(1,nref), [0, 2.9, 0] ];	%[ multiplier ref1, multiplier ref2, offset, linear, quadratic ]
		c = 0:1023;										%linear factor
		cc= c.^2*1e-5;									%quadratic factor is in the range of 1e-5, where all others are order 1
		ccc= [1*ones(size(c)); c; cc]';					%const, lin, qua
	elseif order == 1
		p = [multiplier*ones(1,nref), [-45, 3.08] ];	
		ccc = [1*ones(1,1024); 0:1023]';				%const, lin
	else
		disp( 'order not supported' )
	end
	
	%define the function to minimize
	objective = @(p) 		sum( ( r*p(1:nref)(:) - interp1( ccc*p(1+nref:1+nref+order)(:), s, e, 'linear', 1e10 )'.*wc ).^2 ) - 1e6*sum( p([1:nref,nref+2]).*double( p([1:nref,nref+2]) < 0 ) );
	%options(1)=1 -> verbose, option 5 and 6 -> general Nelder-Mead, option 10 -> max No iterations.
	[p_opt] = fmins( objective, p, [0, 1e-6, 0, 0, 0, 0, 0, 0, 0, 1000] );

	stab_s 	= interp1( ccc*p_opt(nref+1:nref+1+order)(:), s, e, 'linear', 0 );
	p 		= ( p_opt(nref+1:nref+1+order).*[1, 1, 1e-5](1:1+order))(end:-1:1);
	act_con = p_opt(1:nref);
	error	= objective( p_opt );
	if verbose
		figure;
		semilogy( s'.*w', 'b' );
		hold on;
		semilogy( references.*w', 'r' );
		semilogy( references.*w'*act_con(:), 'k' )
		semilogy( stab_s'.*w', 'c' )
		line( [LO, LO], [0.1, 10000], 'color', 'k' )
		line( [HI, HI], [0.1, 10000], 'color', 'k' )
		legend( {"spectrum","references", "stabilized and rebinned", ["reconstruction, ", num2str(act_con(:)')] } );
		disp( objective(p_opt) )
		drawnow;
	end
end
%-----------------------------------------------------------------------------------------------------------------------
%perform FSA using a number of reference spectra and brute-force method
%-----------------------------------------------------------------------------------------------------------------------
function [stab_s, p, act_con, error] = FSA( spectrum, references, order, verbose )
	s = spectrum;
	r = references;
	w = logspace( -2, 0, 1024 );
	w = reshape( w.*double( ( 1:1024 > 100 ) & ( 1:1024 < 900 ) ), size( s ) );
	perform = [0,0,0,0,0,1e9,1];
	Jbest = 1e7;
	%brute force
	for offset = 5:0.1:15
		for linear = 2.85:0.01:2.95 
			for quadratic = 0:5e-6:2.5e-5
				stab_s  = rebin( polyval( [quadratic, linear, offset], 1:1024 ), s(:)' )';
				act_con = bsxfun( @times, r, w ) \ (stab_s.*w);
				error   = sum( ( bsxfun( @times, r, w )*act_con - (stab_s.*w) ).^2 ) + 1e6*any( act_con < 0 );
				if error < Jbest
					Jbest = error;
					perform = [offset, linear, quadratic, act_con(:)', error, any( act_con < 0 )];
				end
				%disp( [offset, linear, quadratic, act_con(:)', error] );
				if 0
					semilogy( w )
					semilogy( s.*w, 'b' );
					hold on;
					semilogy( bsxfun( @times, r, w ), 'r' );
					semilogy( bsxfun( @times, r, w )*act_con, 'k' )
					semilogy( abs( bsxfun( @times, bsxfun( @times, r, w ), act_con' ) ), 'color', 0.6*[1,1,1] );
					legend( {"spectrum","references", "stabilized and rebinned", "reconstruction"} );
				end
			end
		end
	end
	stab_s  = rebin( polyval( perform(3:-1:1) , 1:1024 ), s(:)' )';
	act_con = bsxfun( @times, r, w ) \ (stab_s.*w);
	p       = perform(3:-1:1);
	error   = Jbest;
	if verbose
		disp( perform )
		semilogy( s, 'b' );
		hold on;
		semilogy( r, 'r' );
		semilogy( r*act_con, 'k' )
		legend( {"spectrum","references", "stabilized and rebinned", ["reconstruction, ", num2str(act_con(:)')] } );
		drawnow;
		figure;
		plot( ( bsxfun( @times, r, w )*act_con - (stab_s.*w) ) )
	end
end
%-----------------------------------------------------------------------------------------------------------------------
%find peaks in the spectrum and accurately estimate their location
%-----------------------------------------------------------------------------------------------------------------------
function [stab_s, p, act_con, error] = PEAK_find( spectrum, order, truepeaks, exp_loc, verbose )
	if spectrum > 1025	
		%then its 2048 and some factors need to change some
		[loc2, height2] = findfive( spectrum', 25, 0, truepeaks, verbose, exp_loc );
	else
		[loc2, height2] = findfive( spectrum', 5, 0, truepeaks, verbose, exp_loc );
	end
	[p, s]		   = polyfit( loc2, truepeaks, order );		% calculate scaling factors
	e 	   		   = polyval( p, 1:length( spectrum ) );	% calculate scaled energy axis
	stab_s    	   = max( rebin( e, spectrum' ), 0 );		% rebin to evenly spaced energy axis; prevent negative counts
	act_con 	   = [];									% does not yield activity concentrations
	error		   = sqrt( sum( ( truepeaks - polyval( p, truepeaks ) ).^2 ) );
end
%-----------------------------------------------------------------------------------------------------------------------
%find peaks in the spectrum and accurately estimate their location
%-----------------------------------------------------------------------------------------------------------------------
function [loc2, height2] = findfive( spectrum, N, K, truepeaks, verbose, exp_loc )
	%find the 5 highest peaks in an N-point smoothed logarithmically scaled spectrum:
	%N = 3;
	if isempty( exp_loc )
		%only the length(truepeaks) highest peaks
		%[height, loc] = findpeaks( log10( filtfilt( 1/N*ones(N,1), 1, spectrum ) + 1 )(N:end)', "MinPeakDistance", floor( length( spectrum )/80 ) );
		[height, loc] = findpeaks( filtfilt( 1/N*ones(N,1), 1, log10( spectrum + 1 ) )(N:end-N)' + 1e-15, "MinPeakDistance", floor( length( spectrum )/80 ) );
		%semilogy( loc, 10.^height, '+m' );
		[height, idx] = sort( height, 'descend' );
		loc = N-1 + loc(idx(1:length(truepeaks)));
		height = height(idx(1:length(truepeaks)));
		%sort in ascending order
		[loc, idx] = sort( loc );
		height = height(idx ); 
	else
		%those peaks that are nearest to an expected location
		w = logspace( -2, 0, length(spectrum ) );
		spec2 = ( spectrum.*w );
		%a = log10( filtfilt( 1/N*ones(N,1), 1, spectrum ) + 1 );	%detrend the filtered log-spectrum to make peaks more prominent
		a = filtfilt( 1/N*ones(N,1), 1, spectrum.*w );
		spectrum =  ( a );
		range = floor( length( spectrum )/80 );	%like the "minpeakdistance" parameter, it scales with number of channels
		for i = 1:length( exp_loc )
			[height2, loc2 ] = max( spectrum(exp_loc(i) - range : exp_loc(i) + range) );
			loc(i) = exp_loc(i) - range - 1 + loc2;
			height(i) = height2;
		end
	end
	if verbose
		plot( spec2 ); hold on;
		plot( spectrum );
		plot( loc, height, 'om' );
	end
	%Refine the peak height and location by quadratic interpolation in peak-res to peak+res, where res is energy-dependent:
	seven = length( spectrum ) / 146;	%is close to 7 for 1024 channel spectra, but it scales to 14 for 2048 channel spectra; alternative: seven = 10;
	res = seven*sqrt( truepeaks./662 );	%relative resolution per frequency, assuming 7% at 662 keV 
	res = ceil( res );		    		%results in ~ 7 channels for 511 keV to 14 channels for 2560 keV (for 1024 channel spectra)

	for i=1:length(truepeaks)
		x 			= ( loc(i)-res(i):min( loc(i)+res(i), length(spectrum) ) );	%ROI
		%x 			= x( spectrum(x) > max(spectrum([x(1),x(end)])));			%force 'even' tails to the peak to fit data on 
		xs			= x(1):0.01:x(end);
		p 			= polyfit( x, spectrum(x), 3 );
		[d1, d2] 	= max( polyval( p, xs ) );
		height2(i)  = d1;
		loc2(i)		= xs(d2);
	    if verbose | iscomplex( loc2(i) )
			plot( xs, polyval( p, xs ), 'r' ); 
			plot( loc2(i), interp1( spectrum, loc2(i) ), '+r' ); 
		end
	end	
	if verbose
		drawnow;
		pause;
		clf;
	end
end
%-----------------------------------------------------------------------------------------------------------------------
%calculate the energy resolution
%-----------------------------------------------------------------------------------------------------------------------
function resolution = calcresolution( energy, spectrum, expres )
	for i=1:length(truepeaks)
		res(i)		= 0.10 * sqrt( truepeaks(i)/662 ); 							%expected resolution in channels; slightly overestimated
		x 			= ( loc(i)-res(i):min( loc(i)+res(i), length(spectrum) ) );	%ROI
		xs			= x(1):0.01:x(end);											%finer scale
		y			= interp1( spectrum, xs );									%interpolate spectrum over ROI -> could be band limited interpolation?
		xs 			= xs( y(xs) > max(y([xs(1),xs(end)])));						%'even' tails to the peak to fit data on, reduces ROI width
		p 			= polyfit( xs, y(xs), 2 );
		sigma		= sqrt( 0.5*1/abs(p(1)) );									%in logspace, a peak that was gaussian shaped, becomes quadratic shaped
		FWHM		= 2.35482*sigma;											%gaussian sigma to FWHM relationship
		resolution 	= FWHM / loc(i);											%resolution at energy = truepeaks(i)
		resolution 	= resolution * sqrt( truepeaks(i)/662 );					%correction for not measuring the resolution at 662 keV
	end		
end
%-----------------------------------------------------------------------------------------------------------------------
% display results
%-----------------------------------------------------------------------------------------------------------------------
function show_results( t, data, calibrated_spectra, stabparam, method, LIMS )
	N 		  = size( data, 1 );
	K 		  = size( data, 2 ); 
	LO		  = 100;
	HI		  = 900;
	if method == 1	%images and stabilization parameters
		subplot( 3,1,1 );
		imagesc( bsxfun( @times, data, logspace(-3,0,size(data,1))' ), LIMITS = LIMS(1:2) );
		subplot( 3,1,2 );
		imagesc( bsxfun( @times, calibrated_spectra', logspace(-3,0,size(calibrated_spectra,2))' ), LIMITS = LIMS(1:2) );
		subplot( 3,1,3 );
		plot( stabparam );
		axis( [0, K, LIMS(3:4)] );
	end

	if method == 2	%detrended logarithmic images and stabilization parameters
		figure;
		subplot( 3,1,1 );
		imagesc( ( log10( data./repmat( sum( data ), N, 1 ) + 0.0001 )), LIMITS = LIMS(1:2)  );
		subplot( 3,1,2 );	
		imagesc( ( log10( calibrated_spectra'./repmat( sum( calibrated_spectra' ) , N, 1 ) + 0.0001 )), LIMITS = LIMS(1:2) ); %normalization
		subplot( 3,1,3 );	
		plot( stabparam  );
		axis( [0, K, LIMS(3:4)] );			 
	end		

	if method == 3	%detrended logarithmic figures
		figure;
		plot( ( log10( data./repmat( sum( data(LO:HI,:) ), N, 1 ) + 0.0001 )  ), 'k' ); hold on;
		title( 'raw' )
		%axis( [0,1024,LIMS(1:2)] ) 
		%figure;
		plot( ( log10( calibrated_spectra'./repmat( sum( calibrated_spectra(:,LO:HI)' ) , N, 1 ) + 0.0001 )  ), 'r' );
		title( 'raw (black), stabilized (red)' )
		%axis( [0,1024,LIMS(1:2)] )
	end
end
%-----------------------------------------------------------------------------------------------------------------------
% perform stabilization over a dataset
%-----------------------------------------------------------------------------------------------------------------------
function [calibrated_spectra, stabparam, error, activities] = stabilize( data, method, order, verbose, references, exp_loc, truepeaks )
	N 		  = size( data, 1 ); 		%number of channels
	K 		  = size( data, 2 ); 		%number of spectra
	R 		  = size( references, 2 );	%number of reference spectra
	
	calibrated_spectra 	= zeros( K, N );
	stabparam 			= zeros( K, order+1 );	
	activities 			= zeros( K, R );
		
	for i=1:K
		if sum( data(:,i) ) > 1e4
			if strcmpi( method, 'FSA' )
					[newspectrum, p, act_con, errorloc]  = FSA( data(:,i)(:), references, order, verbose );
				elseif strcmpi( method, 'FSA_optimizers' )
					[newspectrum, p, act_con, errorloc]  = FSA_optimizers( data(:,i)(:), references, order, verbose);
				elseif strcmpi( method, 'PEAK_find' )
					[newspectrum, p, act_con, errorloc]  = PEAK_find( data(:,i)(:), order, truepeaks, exp_loc, verbose );
					act_con = zeros( 1, R );
				else
					disp( 'unrecognized method' ); return
			end
			%show progress
			disp( [i, p, errorloc] ); fflush( stdout ); 
			%save data
			calibrated_spectra( i, : ) = single( newspectrum ); 				
			stabparam( i, : )          = single( p(end:-1:1) );
			activities( i, : )         = single( act_con ); 
			error( i ) 				   = errorloc;
		else
			calibrated_spectra( i, : ) = round( rand(N,1) ); 	
		end
	end

end
%-----------------------------------------------------------------------------------------------------------------------
% read datasets; highly specific for the type of files i worked on...
%-----------------------------------------------------------------------------------------------------------------------
function [data, t, spc] = ReadCSV( filenames )
	% datashape = [N, 1024] or [N, 2048] 	integers or floats
	% tshape = [N, 1]						floats
	% spc = [N,k]							strings
	%
	%file naming convention(s) affect how the data is stored:
	Cobalt   = index( filenames, 'Cobalt' );
	Sodium   = index( filenames, 'Sodium' );
	Vanadium = index( filenames, 'Vanadium' );
	NOAVG    =~index( filenames, 'dailyaverages' ) & index( filenames, 'BG' );
	BG       = index( filenames, 'dailyaverages' );
	OLD      = index( filenames, 'OLD' );
	if Cobalt | Sodium | Vanadium | NOAVG
		fid = fopen( filenames, 'r' );				% read the header lines as text before reading spectrum data 
		spcfile = strsplit( fgetl( fid ), ';' );	% 
		fgetl( fid );								% "real time"
		spcrt = strsplit( fgetl( fid ), ';' );		% values
		fgetl( fid ); 								% "live time"
		spclt = strsplit( fgetl( fid ), ';' ); 		% values
		fgetl( fid ); 								% "dead time"
		spcdt = strsplit( fgetl( fid ), ';' );		% values
		fgetl( fid ); 								% "start date"
		spcdate = strsplit( fgetl( fid ), ';' ); 	% values
		fgetl( fid ); 								% "start time"
		spctime = strsplit( fgetl( fid ), ';' ); 	% values
		fgetl( fid ); 								% "calib coeff ..."
		spcerr = strsplit( fgetl( fid ), ';' ); 	% ERROR msgs on row 13
		fclose( fid );
		%read spectra
		data 	  = dlmread( filenames, ';', 14, 0 );	
		%only non-empty columns:
		idx = ~cellfun( @isempty, spcdate );
		%remove them:
		spcfile = spcfile(idx); 
		spcdate = spcdate(idx);
		spctime = spctime(idx);
		data	= data(:,idx);
		%if error reported in row 13
		idx = ~cellfun( @isempty, spcerr );
		%make those spectra 0
		data(:,idx) = round( rand( size(data,1), nnz(idx) ) ) ;
		%datetime-strings to numbers
		num = datenum( [char(spcdate), char(spctime)], 'dd.mm.yyyyHH:MM:SS' );
		t = num - datenum( '01.01.201500:00:00', 'dd.mm.yyyyHH:MM:SS' );
		spc = strvcat(	cstrcat( [char(spcfile), repmat( ';', size(spcfile,2),1 )] )'(:)', ...
						cstrcat( [char(spcdate), repmat( ';', size(spcdate,2),1 )] )'(:)', ...
						cstrcat( [char(spctime), repmat( ';', size(spctime,2),1 )] )'(:)'  );
	elseif OLD
		fid = fopen( filenames, 'r' );				% read the header lines as text before reading spectrum data 
		spcfile = strsplit( fgetl( fid ), ';' );	% KB...
		fgetl( fid );								% "$DATE_MEA:"
		spcdatetime = strsplit( fgetl( fid ), ';' );% values
		fgetl( fid ); 								% "meas time"
		spclt = strsplit( fgetl( fid ), ';' ); 		% values
		fgetl( fid ); 								% "mca cal"
		spcdt = strsplit( fgetl( fid ), ';' );		% values
		fgetl( fid ); 								% "data"
		spcdate = strsplit( fgetl( fid ), ';' ); 	% values
		fclose( fid );
		%read spectra
		data 	  = dlmread( filenames, ';', 9, 0 );	
		%only non-empty columns:
		idx = ~cellfun( @isempty, spcdate );
		%remove them:
		spcfile 	= spcfile(idx); 
		spcdatetime = spcdatetime(idx);
		data		= data(:,idx);
		%if error reported in row 13
		%idx = ~cellfun( @isempty, spcerr );
		%make those spectra 0
		%data(:,idx) = round( rand( size(data,1), nnz(idx) ) ) ;
		%datetime-strings to numbers
		num = datenum( [char(spcdatetime)], 'dd-mm-yyyy HH:MM' );
		t = num - datenum( '01.01.201300:00:00', 'dd.mm.yyyyHH:MM:SS' );
		spc = strvcat(	cstrcat( [char(spcfile), repmat( ';', size(spcfile,2),1 )] )'(:)', ...
						cstrcat( [char(spcdatetime), repmat( ';', size(spcdatetime,2),1 )] )'(:)'  );
	elseif BG
		fid 	= fopen( filenames, 'r' );			% read the header lines as text before reading spectrum data 
		spcdate = strsplit( fgetl( fid ), ';' );	% date 
		t 		= strsplit( fgetl( fid ), ';' );	% day of year 2015
		fclose( fid );
		data 	  = dlmread( filenames, ';', 2, 0 );
	else 
		disp("Did not recognize file by its name, using fallback routine");
		fid 	= fopen( filenames, 'r' );			% read the header lines as text before reading spectrum data 
		spcdate = strsplit( fgetl( fid ), ';' );	% date 
		t 		= strsplit( fgetl( fid ), ';' );	% day of year 2015
		fclose( fid );
		data 	  = dlmread( filenames, ';', 2, 0 );
	end
end
%-----------------------------------------------------------------------------------------------------------------------
% make a movie of spectrum data
%-----------------------------------------------------------------------------------------------------------------------
function CreateMovie( data )
	%for the movie
	addpath( "C:\\Program Files\\gs\\gs9.14\\bin" )
	system( "PATH=%PATH%;C:\\Program Files\\gs\\gs9.14\\bin" )
	chdir( "C:\\Users\\JanStegenga\\Documents\\Analysis\\movie" )
	%clear directory of .png?
	% system("ffmpeg -r 4 -i frame%04d.png -y -t 00:00:30 movie_rob.avi" ); return %uncomment to make the movie		
	for i = 1:size( data, 2 )
		%create 10 lines with decreasing grayscale intensity
		for j = max([i-9,1]):i
			plot( data(:,j), 'color', [1,1,1] - [1,1,1]*(exp((j-i)/3)) ); hold on;
		end
		title( num2str(i) );
		axis( [0 1024, 0, 12000] );
		drawnow;
		pause(0.1);
		fname = sprintf('frame%04d.png',i); 
		print('-dpng','-r300',fname);
		clf;
	end
end

%-----------------------------------------------------------------------------------------------------------------------
%-----------------------------------------------------------------------------------------------------------------------
% Scripting part
%-----------------------------------------------------------------------------------------------------------------------
%-----------------------------------------------------------------------------------------------------------------------

%true peaks: the energies in keV where peaks should be, run in verbose mode once to check accordance in graphs
%Some prime energies:
%	keV			nuclide
%	510.77, 	208-Tl
%	609.3, 		Bi-214 
%	911.2, 		Ac-228 
%	1460.8, 	K-40 
%	1764.5, 	Bi-214 
%	2614.5, 	208-Tl 
%
truepeaksNa 	= [511 1022 1275 1786 2296];			%Na22
truepeaksCo	    = [811 1173 1333 2506];					%Cobalt 
truepeaksBG		= [510.77 609.3 911.2 1460.8 2614.5];	%background 5x3"
truepeaksBGwell = [510.77, 1460.8, 1764.5, 2614.5];		%background well 4x4"	

%parameters, data files and output files. 
%	exp_loc gives the channel/bin at which to find a truepeak (must have same length). 
%	run in verbose mode to check accordance. leave empty to run without expected locations.
basedir =  "C:\\Users\\"
filenames =  [basedir, "BG.csv"]; 		truepeaks = truepeaksBGwell;	exp_loc = [173 478 572 836];

%set to 1 to run in verbose mode
verbose   = 0; 
%1 = linear stabilization (requires >= 2 peaks), 2 = quadratic stabilization (requires >=3 peaks)
order	  = 2;	
%select one of these
method    = ['FSA', 'FSA_optimizers', 'PEAK_find'](3)

%load data
[data, t, spc] = ReadCSV( filenames );

%load reference spectra from file:
if method[1:3] == 'FSA'
	references = load( [basedir, "references.mat"], "data" );
	references = references.data';
else
	references = []

%stabilize data
[calibrated_spectra, stabparam, error, activities] = stabilize( data, method, order, verbose, references, exp_loc, truepeaks );

%save stabilized data
fileout = [filenames(1:end-4), '-stabilized-', method, '-', num2str( order), '.csv'];
fid 	= fopen( fileout, 'w' ); 
%copy a piece of the 'header' information
for i=1:size( spc,1), fputs( fid, [spc(i,:), "\n"] ); end, fputs( fid, "\n"); fclose( fileout ); 
%write stablization parameters and stabilized spectra
dlmwrite( fileout, stabparam', ';', 0, 0, '-append' );
dlmwrite( fileout, calibrated_spectra', ';', 0, 0, '-append' );

%visualization(s)
show_results( t, data, calibrated_spectra, stabparam, 1, [0 Inf -100 5] );		%do not normalize counts in images
show_results( t, data, calibrated_spectra, stabparam, 2, [-4 -2 -100 5] );		%normalize counts in images

return
