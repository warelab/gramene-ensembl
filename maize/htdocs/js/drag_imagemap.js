function select_stop(evt) { debug_print( 'select_stop' );
  evt = (evt) ? evt : ((event)?event:null)
  document.getElementsByTagName('body')[0].onmousemove = null;
  document.getElementsByTagName('body')[0].onmouseup   = null;

  F = document.forms['panel_form'];
  el = ego( 'other_l' )
  if( el && egH( el )+egW( el ) > 0 ) { /* this is our current rubber band */
    sx = egX( el );
    sy = egY(ego('other_t'));
    hide( el );
    hide(ego('other_t'));
    hide(ego('other_r'));
    hide(ego('other_b'));
    rs( el, 0 , 0 );
    rs( ego('other_t'), 0 , 0 );
    rs( ego('other_r'), 0 , 0 );
    rs( ego('other_b'), 0 , 0 );
    e_x = enx   = egeX(evt);
    e_y = eny   = egeY(evt);
    O    = egi(dragging_id+'_i');
    tl_x = egX(O);
    tl_y = egY(O);
    br_x = tl_x + egW(O)  -1;
    br_y = tl_y + egH(O) -1;
    if( enx < tl_x ) { enx = tl_x; }
    if( eny < tl_y ) { eny = tl_y; }
    if( enx > br_x ) { enx = br_x; }
    if( eny > br_y ) { eny = br_y; }
    W = enx-sx;
    H = eny-sy;
    if( W*W+H*H > N*N-1) {
//  ac(dragging_object,ego('real_l'));
//  ac(dragging_object,ego('real_r'));
//  ac(dragging_object,ego('real_b'));
//  ac(dragging_object,ego('real_t'));
      show(ego('real_l'));
      show(ego('real_t'));
      show(ego('real_r'));
      show(ego('real_b'));

      m2(ego('real_l'), sx,    sy, 1,   H+1 ); m2(ego('real_t'), sx, sy,    W+1, 1   );
      m2(ego('real_r'), enx,   sy, 1,   H+1 ); m2(ego('real_b'), sx, eny,   W+1, 1   );
      chr = F.chr.value
/* Centre point of object */
      ocp = Math.floor( 0.5 * F.elements[dragging_id+'_bp_end'].value + 0.5 * F.elements[dragging_id+'_bp_start'].value );
      ow  = Math.floor( 1.0 * F.elements[dragging_id+'_bp_end'].value - 1.0 * F.elements[dragging_id+'_bp_start'].value + 1 );
      ns  = 1.0 * p2b( dragging_id,  sx-tl_x );
      ne  = 1.0 * p2b( dragging_id, enx-tl_x );
      cp  = Math.floor( 0.5 * ns + 0.5 * ne );
      w   = ne - ns;
      aw  = F.elements['main_width'].value;
      panel_type = F.elements[dragging_id+'_flag'].value;
      switch( panel_type ) {
        case 'bp':
          lnks = new Array(
            new Array( 'View this region in basepair view', cv_URL( { c: cp, zw: w } ) ),
            new Array( 'Centre on this region',             cv_URL( { c: cp }       ) )
          );
          break;
        case 'cv':
          lnks = new Array(
//            new Array( 'View this region in basepair view', cv_URL( { c: cp, zw: w } ) ),
            new Array( 'Zoom into this region on detailed view', cv_URL( { c: cp, w: w  } ) ),
            new Array( 'Centre on this region',             cv_URL( { c: cp        } ) )
          );
      }
      var S = seq_region( dragging_id ) +': '+ns+ " - " + ne;
      show_zmenu( S, e_x, e_y, lnks, 0 );
      return void(0);
    } else {
// This is a click 
      image_map = ego( dragging_id + '_i_map' );
// If we have an image map lets see if it is in a feature
/*
      if(evt.ctrlKey) {
        chr = F.chr.value
        cp  = p2b( dragging_id, enx-tl_x )
        aw  = F.elements['main_width'].value;
        // This is the centering action!! 
        URL = '';
        panel_type = F.elements[dragging_id+'_flag'].value;
        switch( panel_type ) {
          case 'bp': URL = cv_URL( { c: cp } ); break;
          case 'cv': URL = cv_URL( { c: cp } ); break;
        }
        alert( URL )
        return void(0);
      } else if(evt.shiftKey) {
        // Make this a zoom in action.... 
        chr = F.chr.value 
        cp  = p2b( dragging_id, enx-tl_x )
        aw  = F.elements['main_width'].value;
        URL = '';
        panel_type = F.elements[dragging_id+'_flag'].value;
        switch( panel_type ) {
          case 'bp': URL = cv_URL( { c: cp, zw: F.elements['bp_width'].value/2   } ); break;
          case 'cv': URL = cv_URL( { c: cp,  w: F.elements['main_width'].value/2 } ); break;
        }
        alert( URL )
        return void(0);
      }
*/
      if(image_map && !(evt.altKey)) {
        L = image_map.areas.length;
        for(i=0;i<L;i++) {
          zmn = 'zmenu_'+dragging_id+'_'+i
          A = image_map.areas[i];
          pts = A.coords.split(/\D+/);
          if( A.shape=='poly' ? in_poly( pts, enx-tl_x, eny-tl_y ) : ( A.shape=='circle' ? in_circle( pts, enx-tl_x, eny-tl_y ) : in_rect( pts, enx-tl_x, eny-tl_y ) ) ) {
            if(A.onclick) {
              CLICK_X  = e_x;
              CLICK_Y  = e_y;
              ZMENU_ID = zmn;
              MOUSE_UP = 1;
              A.onclick();
              MOUSE_UP = 0;
            } else if(A.onmouseover) {
              CLICK_X  = e_x;
              CLICK_Y  = e_y;
              ZMENU_ID = zmn;
              MOUSE_UP = 1;
              A.onmouseover();
              MOUSE_UP = 0;
            } else {
              if( A.title.substr(  0, 6 ) == 'About:' ) {
                lnks = new Array(
                  new Array( "Genes were annotated by the Ensembl automatic analysis pipeline using either a GeneWise model from a human/vertebrate protein, a set of aligned human cDNAs followed by GenomeWise for ORF prediction or from Genscan exons   supported by protein, cDNA and EST evidence. GeneWise models are further combined with available aligned cDNAs to annotate UTRs." ),
                  new Array('Further help',A.href)
                );
              } else {
/*
                if( A.href && A.href != 'javascript:void(0)') {
                  lnks = new Array(
                    new Array( 'test' ),
                    new Array('Further information...',A.href)
                  );
                } else {
                  lnks = new Array(
                    new Array( 'test' )
                  );
                }
*/
                if( A.href ) {
                  location.href = A.href; return;
                }
              }
              show_zmenu( A.title, e_x, e_y, lnks, zmn );
            }
            return void(0);
          }
        }
      }
// If not give the appropriate centering dialog...
      chr = F.chr.value
      cp  = p2b( dragging_id, enx-tl_x )
      aw  = 1.0*F.elements['main_width'].value; 
      zw  = 1.0*F.elements['bp_width'].value; 
      panel_type = F.elements[dragging_id+'_flag'].value;
      switch( panel_type ) {
        case 'bp': 
          lnks = new Array(
//        new Array( 'Centre basepair view here', '/' ),
            new Array( 'Zoom in x5',    cv_URL({ c: cp, zw: zw/5 } ) ),
            new Array( 'Zoom in x2',    cv_URL({ c: cp, zw: zw/2 } ) ),
            new Array( 'Centre',        cv_URL({ c: cp, zw: zw } ) ),
            new Array( 'Zoom out x2',   cv_URL({ c: cp, zw: zw*2 } ) ),
            new Array( 'Zoom out x5',   cv_URL({ c: cp, zw: zw*5 } ) )
          );
          break;
        case 'cv':
          lnks = new Array(
//        new Array( 'Centre basepair view here', '/' ),
            new Array( 'Zoom in x10',   cv_URL({ c: cp, w: aw/10 } ) ),
            new Array( 'Zoom in x5',    cv_URL({ c: cp, w: aw/5 }  ) ),
            new Array( 'Zoom in x2',    cv_URL({ c: cp, w: aw/2 }  ) ),
            new Array( 'Centre',        cv_URL({ c: cp          }  ) ),
            new Array( 'Zoom out x2',   cv_URL({ c: cp, w: aw*2 }  ) ),
            new Array( 'Zoom out x5',   cv_URL({ c: cp, w: aw*5 }  ) ),
            new Array( 'Zoom out x10',  cv_URL({ c: cp, w: aw*10 } ) )
          );
      }
      show_zmenu( seq_region( dragging_id ) +': '+ cp, e_x, e_y, lnks );
    }
  }
  return void(0);
}
