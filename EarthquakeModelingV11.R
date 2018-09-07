##########################################################################################
###
### EarthquakeModeling.R
### MATH3800: Stochastic Processes
###
### Final
### Project:  Earth Quake Prediction Modeling
###
### Authors:  Matt Bennett (bennetts4@wit.edu)
###           Kai Yuen Fong
###           Simranjit Singh
###           Ichiryu Nakashima
###
### Summary:  Simulation of major earthquake occurances in the Southern California region
###           using data files obtained from CalTech: Division of Geological and Planetary
###           Sciences. (Simulation of other data sets may require changes to the parsing
###           and structuring of the global data variable.)
###           
##########################################################################################

### Required for map display functions:
###   Tools -> Install Packages => ggmap from repository
###   Tools -> Install Packages => mapproj from repository
###
### NOTE: LIBRARY AND PROVIDED FUNCTIONS USE LON, LAT
library( ggmap )
library( mapproj )

### NOTE: DATA FILES MUST BE IN YOUR WORKING DIRECTORY
workingDirectory = "C:/Users/bennetts4/Google Drive/zCompleted_Courses/MATH3800_StochasticProcesses/Stochastics_Project/DataFiles";
setwd(workingDirectory);
# Number of years to simulate
t = 50;

# Boundaries of data set
startYear = 1932;    endYear = 2016;
minLat = 31.999;     maxLat = 37.001;
minLon = -122.001;   maxLon = -113.999;

# Information on how to structure the data
minMag = 3;
mainShockPercent = 0.95;
rows = 8;      columns = 8;
subrows = 4;   subcolumns = 4;

# Conversions of structural data
startDate = as.Date( paste( startYear, "-01-01", sep="" ) );
endDate = as.Date( paste( endYear, "-01-01", sep="" ) );
deltaLat = ( maxLat - minLat ) / columns;
deltaLon = ( maxLon - minLon ) / rows;
subDeltaLat = deltaLat / subcolumns;
subDeltaLon = deltaLon / subrows;

### Parse earthquake data from text files startYear.txt through endYear.txt in 
### working directory
###
### Global      rawData      data frame      
parseData = function(){
    data = c();
    for ( year in startYear:endYear ){
        fileData = read.table( paste(year, ".txt", sep = ""), header = TRUE );
        data = rbind( data, fileData );
    }
    rawData <<- data; # make rawData available as global variable in R
    print( "Data parsed from files." );
    return( rawData );
}

### Clean and structure data for use with R functions.
###
### Global      data        data frame      provided data with unused data and
###                                         information removed, data types altered,
###                                         and information for calculations added
###                                         in new columns
structureData = function( data ){
    data = subset( data, ( data$MAG > minMag ) );
    data$YYY.MM.DD = as.Date( data$YYY.MM.DD );
    data$DATE = data$YYY.MM.DD;

    data$YYY.MM.DD = NULL;   data$HH.mm.SS.ss = NULL;   data$DEPTH = NULL;
    data$ET = NULL;          data$GT = NULL;            data$NGRM = NULL;
    data$M = NULL;           data$Q = NULL;             data$NPH = NULL;
    
    data = createGridData( data );
    data = createEventData( data );
    data = createQuantileData( data );
    data = createLambdaData( data );
    data <<- data # make structured data available as global variable in R
    print( "Data structured and cleaned." );
    return( data );
}



### Add row and column, and subrow, subcolumn, information to the data set.
###
### Global      data        data frame      provided data with additional
###                                         information for region
createGridData = function( data ){
    data$COLUMN = ceiling( ( data$LAT - minLat ) / deltaLat );
    data$ROW = ceiling( ( maxLon - data$LON ) / deltaLon );
    data$SUBCOLUMN = ceiling( ( data$LAT - minLat ) / subDeltaLat ) %% subcolumns + 1;
    data$SUBROW = ceiling( ( maxLon - data$LON ) / subDeltaLon ) %% subrows + 1;
    return( data );
}

### Access all the occurances in the provided data set of the given row and column 
### or of the row, column, subrow, and subcolumn.
###
### Local       data        data frame      all data from input data within the
###                                         specified section of the overall
###                                         region           
getCellData = function( data, column, row, subcolumn, subrow ){
    data = subset( data, ( data$COLUMN == column ) & ( data$ROW == row ) );
    if ( missing( subcolumn ) | missing( subrow ) ){
        return( data );
    } else{
        data = subset( data, ( data$SUBCOLUMN == subcolumn ) & 
                           ( data$SUBROW == subrow ) );
        return( data );
    }
}

### Determine event defining earthquakes and cluster nearby foreshocks and after
### shocks together with the event using a event idenifier.
### 
### Local       eventData   data frame      all data from input data with event
###                                         column added with subsets of matching 
###                                         interger event identifiers
###                                         representing a clustered event.
###
createEventTriggers = function( data ){
    eventData = c();
    for ( row in 1:rows ){ for ( column in 1:columns ){
        for ( r in 1:subrows ){ for ( c in 1:subcolumns ){
            subcell = getCellData( data, column, row, c, r );
            threshold = quantile( subcell$MAG, prob = mainShockPercent )
            mainShock = subset( subcell, ( subcell$MAG >= threshold ) );
            eventData = rbind( eventData, mainShock );
        }}
    }}
    eventData = eventData[ !( 
            duplicated( eventData$DATE ) & 
            duplicated( eventData$COLUMN ) & duplicated( eventData$ROW ) & 
            duplicated( eventData$SUBCOLUMN ) & duplicated( eventData$SUBROW ) 
    ), ];
    eventData = eventData[ with( eventData, order( DATE, MAG) ), ];
    eventData$EVENT = c( 1:nrow( eventData ) );
    return( eventData );
}

###
###
###
createEventData = function( data ){
    data$EVENT = -1;
    eventData = createEventTriggers( data );
    data = data[ with( data, order( DATE, MAG) ), ];
    for ( e in 1:nrow( eventData ) ){
        data$EVENT[ data$EVID == eventData$EVID[ e ] ] = eventData$EVENT[ e ];
        subcell = getCellData( data, eventData$COLUMN[ e ], eventData$ROW[ e ], 
                               eventData$SUBCOLUMN[ e ], eventData$SUBROW[ e ] );
        if ( nrow( subcell ) != 0 ){
            for ( i in 1:nrow( subcell ) ){
                if ( ( subcell$DATE[ i ] >= eventData$DATE[ e ] - 3 ) & 
                     ( subcell$DATE[ i ] <= eventData$DATE[ e ] + 3 ) ){
                    data$EVENT[ data$EVID == subcell$EVID[ i ] 
                                ] = eventData$EVENT[ e ];
                }   }   }
    }
    return( data );
}

###
###
###
getEventTriggers = function( data, column, row ){
    eventData = c();
    if ( missing( column ) || missing( row ) ){
        for ( c in 1:columns ){ for ( r in 1:rows ){
            cell = getCellData( data, c, r );
            t = getEventTriggers( cell, c, r );
            if ( is.null( eventData ) ){
                eventData = t;
            }else {
                eventData = rbind( eventData, t );
            }
        }}
    }else {
        cell = getCellData( data, column, row );
        u = as.data.frame( unique( cell$EVENT ) );
        u = subset( u, u[ , 1] > 0 );
        for ( i in 1:nrow( u ) ){
            e = subset( cell, cell$EVENT == u[ i, ] );
            if ( nrow( e ) != 0 ){
                t = subset( e, e$MAG == max( e$MAG ) );
                if ( is.null( eventData ) ){
                    eventData = t;
                }else {
                    eventData = rbind( eventData, t );
                }
            }
        }
    }
    return( eventData );
}

###
###
###
createQuantileData = function( data ){
    data$QUANT = 0;
    for ( row in 1:rows ){ for ( column in 1:columns ){
        cell = getCellData( data, column, row );
        if ( nrow( cell ) != 0 ){
            cell$QUANT = quantile( cell$MAG, prob = mainShockPercent );
            for ( i in 1:nrow( cell ) ){
                data$QUANT[ data$EVID == cell$EVID[ i ] ] = cell$QUANT[ i ];
            }
        }
    }}
    return( data );
}

###
###
###
createLambdaData = function( data ){
    data$LAMBDA = 0;
    for ( row in 1:rows ){ for ( column in 1:columns ){
        cell = getCellData( data, column, row );
        if ( nrow( cell ) != 0 ){
            cell$LAMBDA = getCellLambda( cell );
            for ( i in 1:nrow( cell ) ){
                data$LAMBDA[ data$EVID == cell$EVID[ i ] ] = cell$LAMBDA[ i ];
            }
        }
    }}
    return( data );
}

###
###
###
getCellLambda = function( cell ){
    if ( nrow( cell ) == 0 ){
        ilambda = 0;
    } else{
        events = getEventTriggers( cell );
        if ( is.null( events ) ){
            ilambda = 0;
        }else if ( nrow( events ) > 4 ){
            ilambda = 365 / (as.numeric( max( events$DATE ) - min( events$DATE ) ) / nrow( events ) - 1 );
        } else{
            ilambda = 365 / as.numeric( endDate - startDate ) / nrow( events );
        }
    }
    return( ilambda );
}



###
###
###
getCellQMean = function( cell ){
    cell = getEventTriggers( cell );
    cell = subset( cell, cell$MAG > cell$QUANT );
    if ( is.null( cell ) || (nrow( cell ) == 0) ){
        avg = 0;
    } else{
        avg = mean( cell$MAG );
    }
    return( avg );
}

###
###
###
getCellMag = function( cell ){
    cell = getEventTriggers( cell );
    if ( is.null( cell ) || (nrow( cell ) == 0) ){
        mag = NULL;
    }else{
        Q = cell$QUANT[ 1 ];
        EM = getCellQMean( cell );
        if ( (Q == 0) || (EM == 0) ){
            mag = NULL;
        } else{
            R = runif( 1, min = 0, max = 1 );
            mag = Q - ( EM - Q ) * log( 1 - R );
        }
    }
    return( mag );
}

getPois = function( cell ){
    lambda = getCellLambda( cell );
    v = runif( rpois( 1, lambda * t ), min = 0, max = t );
    v = sort( v );
    # plot( stepfun( v, 0:length( v ) ) );
    # print( length( v ) );
    # print( v );
    return( v );
}

###
###
###
plotGridData = function( data ){
    par( mfrow = c( columns, rows ), oma = c( 1, 1, 0, 0 ),
         mar = c( 1, 1, 1, 0 ), tcl = -0.1, mgp = c( 0, 0, 0 ) 
    );
    for ( r in 1:rows ){ for ( c in 1:columns ){
        d = subset( data, ( data$COLUMN == c ) & ( data$ROW == r ) );
        if( nrow( d ) != 0 ){
            p = plot( d$DATE, d$MAG, ylim = c( minMag - 0.2, 8 ), 
                      xlim = c( startDate, endDate ) );
        } else{
            plot( 0, type="l", ylim = c( minMag - 0.2, 8 ), 
                  xlim = c( startDate, endDate ) );
        }
    }}
}

###
###
###
###
###

###
###
###
plotSubGridData = function( data, column, row ){
    par( mfrow = c( subcolumns, subrows ), oma = c( 1, 1, 0, 0 ),
         mar = c( 1, 1, 1, 0 ), tcl = -0.1, mgp = c( 0, 0, 0 ) 
    );
    for ( r in 1:subrows ){ for ( c in 1:subcolumns ){
        s = getCellData( data, column, row, c, r );
        if( nrow( s ) != 0 ){
            p = plot( s$DATE, s$MAG );
        } else{
            plot( 0, type="l" );
        }
    }}
}

###
###
###
plotCellEvents = function( data, column, row ){
    cell = getCellData( data, column, row );
    if ( nrow( cell ) != 0 ){
        triggers = getEventTriggers( cell );
        plot( triggers$DATE, triggers$MAG, 
              xlim = c( startDate, endDate ),
              ylim = c( minMag, 8 ) 
        );
    }else {
        plot( 0, type="l", xlim = c( startDate, endDate ),
              ylim = c( minMag, 8 ) );
    }
}

###
###
###
plotGridEvents = function( data ){
    par( mfrow = c( columns, rows ), oma = c( 1, 1, 0, 0 ),
         mar = c( 1, 1, 1, 0 ), tcl = -0.1, mgp = c( 0, 0, 0 ) 
    );
    for ( r in 1:rows ){ for( c in 1:columns ){ 
        plotCellEvents( data, c, r );
    }}
}

###
###
###
plotPoisMag = function( cell ){
    et = getPois( cell );
    if ( is.null( et ) || (length( et ) == 0) ){
        plot( 0, type="l", xlim= c( 0, t ), ylim = c( 3, 8 ) );
        grid( t/5, 5 );
        return( NULL );
    }
    mag = c();
    for ( e in 1:length( et ) ){
        if ( length( mag ) == 0 ){
            mag = getCellMag( cell );
        } else{
            mag = cbind( mag, getCellMag( cell ) );
        }
    }
    if ( (length( et ) == 0) || is.null( mag ) ){
        plot( 0, type="l", xlim= c( 0, t ), ylim = c( 3, 8 ) );
        grid( t/5, 5 );
        return( NULL );
    } else{
        plot( et, mag, xlim= c( 0, t ), ylim = c( 3, 8 ) );
        grid( t/5, 5 );
    }
    df = c(); 
    if ( length( mag ) == 1 ){
        df$MAG = mag[ 1 ]; df$YEAR = et;
    } else{
        df$MAG = mag[ 1, ]; df$YEAR = et;
    }
    return( df );
}

###
###
###
plotPoisMagGrid = function( data ){
    par( mfrow = c( columns, rows ), oma = c( 1, 1, 0, 0 ),
         mar = c( 1, 1, 1, 0 ), tcl = -0.1, mgp = c( 0, 0, 0 ) 
    );
    simData = c();
    for ( r in 1:rows ){ for( c in 1:columns ){
        cell = getCellData( data, c, r );
        result = c();
        result = plotPoisMag( cell );

        if ( is.null( result ) || (length( result ) == 0) ){
            #pass
        } else if ( is.null( simData ) || (length( simData ) == 0) ){
            result$COLUMN = c;
            result$ROW = r;
            simData = as.data.frame( result );
        } else{
            result$COLUMN = c;
            result$ROW = r;
            simData = merge.data.frame( simData, as.data.frame( result ), all = TRUE );
        }
        simData <<- simData;
    }}
    print( paste( "Simulation of ", t, " years complete.", sep="" ) );
    View( simData );
}

###
###
###
mapGrid = function( simData ){
    mg = c();
    for ( r in 1:rows ){ for ( c in 1:columns ){
        cell = getCellData( simData, c, r );
        m = c();
        if( is.null(cell) || (nrow(cell) == 0 )){
            m$MAG = 0;
            m$FREQ = 0;
        }else{
            m$MAG = max( cell$MAG );
            m$FREQ = (nrow(cell)/nrow(simData));
        }
        m$ROW = r;
        m$COLUMN = c;
        m$maxLAT = (maxLat - (m$ROW-1)*deltaLat);
        m$minLAT = (maxLat - m$ROW*deltaLat);
        m$maxLON = (minLon + m$COLUMN*deltaLon);
        m$minLON = (minLon + (m$COLUMN-1)*deltaLon);
        m = as.data.frame( m );
        if( is.null(mg)||(nrow(mg)==0)){
            mg = m;
        }else{
            mg = rbind( mg, m );
        }
    }}
    print( mg );
    mg <<- mg;
    gm = get_map( "California", zoom = 6, maptype="terrain", source="google" );
    ggmap( gm, base_layer = ggplot( x=minLON, y=minLAT, data=mg ) ) + 
        geom_rect( aes( xmin=minLON, xmax=maxLON, ymin=minLAT, ymax=maxLAT, 
                        fill = I("firebrick1"), alpha=(MAG^2)/64 ) ) + 
        geom_rect( aes( xmin=minLON, xmax=maxLON, ymin=minLAT, ymax=maxLAT, 
                        fill = I("blue"), alpha=(FREQ*2.5) ) );
}

###
###
###
mapData = function( data, zoom, size ){
    if ( missing( zoom ) ){
        zoom=7;
    }
    qmplot( LON, LAT, data=data, size=I(size),
            color=I("red"), zoom=zoom,
            source="google", maptype="terrain" );
}

###
###
###
run = function(){
    setwd( workingDirectory );
    rawData = parseData();
    data = structureData( rawData );
    plotPoisMagGrid( data );
}

###
###
###
quickRun = function(){
    data = structureData( rawData );
    plotPoisMagGrid( data );
}