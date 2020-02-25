//
// Vortexje -- Lifting surface.
//
// Copyright (C) 2012 - 2014 Baayen & Heinz GmbH.
//
// Authors: Jorn Baayen <jorn.baayen@baayen-heinz.com>
//

#include <vortexje/lifting-surface.hpp>
#include <vortexje/parameters.hpp>

using namespace std;
using namespace Eigen;
using namespace Vortexje;

/**
   Constructs an empty LiftingSurface.
*/
LiftingSurface::LiftingSurface(const string &id) : Surface(id) {
    beam_is_there = false;
    beam_location = 0.;
}

/**
   Returns the number of chordwise nodes.

   @returns The number of chordwise nodes.
*/
int 
LiftingSurface::n_chordwise_nodes() const {
    return (int) upper_nodes.rows();
}

/**
   Returns the number of chordwise panels.
   
   @returns The number of chordwise panels.
*/
int 
LiftingSurface::n_chordwise_panels() const {
    return (int) upper_panels.rows()+lower_panels.rows();
}

/**
   Returns the number of spanwise nodes.

   @returns The number of spanwise nodes.
*/
int 
LiftingSurface::n_spanwise_nodes() const {
    return (int) upper_nodes.cols();
}

/**
   Returns the number of spanwise panels.
   
   @returns The number of spanwise panels.
*/
int 
LiftingSurface::n_spanwise_panels() const {
    return (int) upper_panels.cols();
}

/**
   Returns the index'th trailing edge node.
  
   @param[in]   index   Trailing edge node index.
   
   @returns The node number of the index'th trailing edge node.
*/
int 
LiftingSurface::trailing_edge_node(int index) const {
    return upper_nodes(upper_nodes.rows() - 1, index);
}

/**
   Returns the index'th upper trailing edge panel.
  
   @param[in]   index   Trailing edge panel index.
   
   @returns The panel number of the index'th upper trailing edge panel.
*/
int
LiftingSurface::trailing_edge_upper_panel(int index) const {
    return upper_panels(upper_panels.rows() - 1, index);
}

/**
   Returns the index'th lower trailing edge panel.
  
   @param[in]   index   Trailing edge panel index.
   
   @returns The panel number of the index'th lower trailing edge panel.
*/
int
LiftingSurface::trailing_edge_lower_panel(int index) const {
    return lower_panels(lower_panels.rows() - 1, index);
}

/**
   Finishes the set up of the trailing edge.  
   
   This function computes the trailing edge bisectors, the initial wake strip normals, and terminates the neigbor relationships
   along the trailing edge.
*/
void
LiftingSurface::finish_trailing_edge() {
    // Compute trailing edge bisectors and normals to the initial wake strip surface:
    trailing_edge_bisectors.resize(n_spanwise_nodes(), 3);
    wake_normals.resize(n_spanwise_nodes(), 3);
    
    if (n_chordwise_nodes() > 1) {
        for (int i = 0; i < n_spanwise_nodes(); i++) {
            // Compute bisector: 
            Vector3d upper = nodes[upper_nodes(upper_nodes.rows() - 1, i)] - nodes[upper_nodes(upper_nodes.rows() - 2, i)];
            Vector3d lower = nodes[lower_nodes(lower_nodes.rows() - 1, i)] - nodes[lower_nodes(lower_nodes.rows() - 2, i)];
            
            upper.normalize();
            lower.normalize();
            
            Vector3d trailing_edge_bisector = upper + lower;
            trailing_edge_bisector.normalize();
            
            trailing_edge_bisectors.row(i) = trailing_edge_bisector;
            
            // Compute normal to the initial wake strip surface, spanned by the bisector and by the span direction:
            int prev_node, next_node;
            
            if (i > 0)
                prev_node = trailing_edge_node(i - 1);
            else
                prev_node = trailing_edge_node(i);
            
            if (i < n_spanwise_nodes() - 1)
                next_node = trailing_edge_node(i + 1);
            else
                next_node = trailing_edge_node(i);
                
            Vector3d wake_normal(0, 0, 0);
            
            if (prev_node != next_node) {
                Vector3d span_direction = nodes[next_node] - nodes[prev_node];
                
                wake_normal = span_direction.cross(trailing_edge_bisector);
                wake_normal.normalize();
            }
            
            wake_normals.row(i) = wake_normal;
        }
        
    } else {
        // No bisector information available:
        trailing_edge_bisectors.setZero();
        wake_normals.setZero();
    }
    
    // Terminate neighbor relationships across trailing edge.
    for (int i = 0; i < n_spanwise_panels(); i++)
        cut_panels(trailing_edge_upper_panel(i), trailing_edge_lower_panel(i));
}

/**
   Build the nodes for the beam of the lifting surface. The beam is supposed to 
   be located at the quarter of the chord.
*/
void
LiftingSurface::create_beam_nodes() {
    /*
    int node_id  ;
    cout << " !!!!!!!!!!!!!!!!!!! FINISH_BEAM_NODES "<< (int) upper_nodes.rows() << " " << (int)lower_nodes.rows()<< endl;
    for (int i=0; i<(int) upper_nodes.rows(); i++ ){
        node_id = upper_nodes(i,0) ;
        cout << node_id << " " << nodes[node_id][0]<< " " << nodes[node_id][1]<< " " << nodes[node_id][2] << endl;
    } 
    cout << " !!!!!!!!!!!!!!!!!!! \n";
    for (int i=0; i<(int) lower_nodes.rows(); i++ ){
        node_id = lower_nodes(i,0) ;
        cout << node_id << " " << nodes[node_id][0]<< " " << nodes[node_id][1]<< " " << nodes[node_id][2] << endl;
    } */
    
//    assert (this->beam_is_there);
    
    beam_nodes.resize            (n_spanwise_nodes()  );
    beam_nodes_collocation.resize(n_spanwise_nodes()-1);
        
    for (int i = 0; i < n_spanwise_nodes(); i++) {
        int ind_leading   = upper_nodes ( 0                    , i ) ;
        int ind_trailing  = upper_nodes ( upper_nodes.rows()-1 , i ) ;
        
        beam_nodes[i] = nodes[ind_leading] + (nodes[ind_trailing] - nodes[ind_leading]) * beam_location  ;
//        cout << beam_nodes[i](0) << " "<< beam_nodes[i](1) << " "<< beam_nodes[i](2) << endl;
        
        if (i>0) 
            beam_nodes_collocation[i-1] = ( beam_nodes[i] + beam_nodes[i-1] ) * 0.5  ;
    }
/*    
    cout << "----------------------------------------------" <<endl;
    double omega = 2.0 * 7.5010760592212296557026381848104 ; 
    double amplitude = 0.5; // amplitude of oscillation in degrees
    amplitude *= 3.141592653589793238462643383279502884 / 180.0 ;
    double sinus = cos(0.*omega);
    double rotation_value = amplitude*sinus;
    Eigen::Vector3d vectorX = Eigen::Vector3d(1.,0.,0.);
    Eigen::Vector3d vectorP ,locnode;
    Eigen::Transform<double, 3, Eigen::Affine> transformation ;
    double rmin = beam_nodes[0].norm();
    double rmax = beam_nodes[n_spanwise_nodes()-1].norm();
    double rfact,rcurrent,amplitude_local;
    vector_aligned<Eigen::Vector3d> vector_beam;
    vector_beam.resize(n_spanwise_nodes());
    for (int i=0; i<n_spanwise_nodes(); i++ ){
        if (i==0) 
            vector_beam[i] =   beam_nodes[i+1] - beam_nodes[i]   ;
        else if (i==n_spanwise_nodes()-1) 
            vector_beam[i] =   beam_nodes[i]   - beam_nodes[i-1]  ;
        else 
            vector_beam[i] = ( beam_nodes[i+1] - beam_nodes[i-1] ) * 0.5 ;
        
        vector_beam[i].normalize();
        rcurrent = beam_nodes[i].norm();
        rfact    = (rcurrent-rmin)/(rmax-rmin);
        amplitude_local = rfact * rotation_value;
        vectorP = vector_beam[i].cross(vectorX);
        vectorP.normalize();
        transformation = AngleAxisd(amplitude_local, vectorP)  ;
        locnode = transformation * beam_nodes[i];
        cout << vectorP[0] << " " << vectorP[1] << " " << vectorP[2] << "   ---  "<<amplitude_local<< endl;
        cout << beam_nodes[i](0) << " "<< beam_nodes[i](1)<< " "<< beam_nodes[i](2)<< " // "<< locnode(0) << " "<< locnode(1)<< " "<< locnode(2)<< endl;
    }
    cout << "----------------------------------------------" <<endl;*/
}

/**
   Impose an arbitrary motion to the beam of the airfoil and move the airfoil
   accordingly
*/
void
LiftingSurface::motion_beam_nodes(const double &time, const double &dt,
                                  vector_aligned< Eigen::Transform<double, 3, Eigen::Affine> > &transforms_TE) {
    cout << "Surface "<< this->id << ": Deforming the beam of the blade..." <<endl;
    
// Flaping motion of the lifting surface, vector_beam is parallel to the beam    
    vector_aligned<Eigen::Vector3d> vector_beam;
    vector_beam.resize(n_spanwise_nodes());
    
    // Total number of nodes per radial section (-2 because leading and trailing nodes)
    // upper: 0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30   (31)
    // lower: 0 58 57 56 55 54 53 52 51 50 49 48 47 46 45 44 43 42 41 40 39 38 37 36 35 34 33 32 31 30      (30)
//    int numnodessection = (int)upper_nodes.rows() + (int)lower_nodes.rows() - 2 ;
    
    // The frequency of torsion is twice as big as the rotation speed
    double omega = 2.0 * 7.5010760592212296557026381848104 ; 
    double amplitude = 0.5; // amplitude of oscillation in degrees
    amplitude *= 3.141592653589793238462643383279502884 / 180.0 ;
    double sinus = cos(time*omega);
    double rotation_value = amplitude*sinus;
    
    Eigen::Transform<double, 3, Eigen::Affine> transformation ;
    transforms_TE.resize(n_spanwise_nodes());
    
    double rmin = beam_nodes[0].norm();
    double rmax = beam_nodes[n_spanwise_nodes()-1].norm();
    double rfact,rcurrent,amplitude_local;
    Eigen::Vector3d vectorX = Eigen::Vector3d(1.,0.,0.);
    Eigen::Vector3d vectorP ;
        
// Loop on the radial sections, to apply the motion on each section
    for (int i=0; i<n_spanwise_nodes(); i++ ){
        if (i==0) 
            vector_beam[i] =   beam_nodes[i+1] - beam_nodes[i]   ;
        else if (i==n_spanwise_nodes()-1) 
            vector_beam[i] =   beam_nodes[i]   - beam_nodes[i-1]  ;
        else 
            vector_beam[i] = ( beam_nodes[i+1] - beam_nodes[i-1] ) * 0.5 ;
                
        vector_beam[i].normalize();
        vector<int> section_nodes;
        for (int j=0; j<(int) upper_nodes.rows(); j++ )
            section_nodes.push_back ( upper_nodes(j,i) );
        for (int j=1; j<(int) lower_nodes.rows()-1; j++ )
            section_nodes.push_back ( lower_nodes(j,i) );
        
        rcurrent = beam_nodes[i].norm();
        rfact    = (rcurrent-rmin)/(rmax-rmin);
        amplitude_local = rfact * rotation_value;
        vectorP = vector_beam[i].cross(vectorX);
        
//        transformation = Translation3d(beam_nodes[i]) * AngleAxisd(amplitude_local, vector_beam[i]) * Translation3d(-beam_nodes[i]) ;
        transformation = AngleAxisd(amplitude_local, vectorP)  ;
        
//        cout << vector_beam[i][0] << " "<< vector_beam[i][1] << " "<< vector_beam[i][2] << " "<< endl;
        
        for (int j=0; j<(int) section_nodes.size(); j++ ){
            nodes[ section_nodes[ j ] ] = transformation * nodes[ section_nodes[ j ] ];
            
            // Update the transform operation for the trailing edge nodes
            // -> required to update the position for the 1st advected wake nodes
            // see body -> set_beam_motion_lifting_surfaces()
            if ( section_nodes[j] == trailing_edge_node(i) )
                transforms_TE[i] = transformation ;
        }
        Vector3d trailing_edge_bisector = trailing_edge_bisectors.row(i);
        trailing_edge_bisectors.row(i) = transformation.linear() * trailing_edge_bisector;
        
        Vector3d wake_normal = wake_normals.row(i);
        wake_normals.row(i) = transformation.linear() * wake_normal;
    }
//    cout << "-------------------------"<<endl;
    
// Update the geometry (normals, panel collocation points, etc.)    
    compute_geometry();
    if (this->beam_is_there)
        create_beam_nodes();
}

/**
   Transforms this lifting surface.
   
   @param[in]   transformation   Affine transformation.
*/
void LiftingSurface::transform(const Eigen::Transform<double, 3, Eigen::Affine> &transformation) {
    // Call super:
    this->Surface::transform(transformation);
    
    // Transform bisectors and wake normals:
    transform_bisectors_wake_normals(transformation);
}

void LiftingSurface::transform_bisectors_wake_normals(const Eigen::Transform<double, 3, Eigen::Affine> &transformation) {
    // Transform bisectors and wake normals:
    for (int i = 0; i < n_spanwise_nodes(); i++) {
        Vector3d trailing_edge_bisector = trailing_edge_bisectors.row(i);
        trailing_edge_bisectors.row(i) = transformation.linear() * trailing_edge_bisector;
        
        Vector3d wake_normal = wake_normals.row(i);
        wake_normals.row(i) = transformation.linear() * wake_normal;
    }
}

/**
   Returns the wake emission velocity, at the node_index'th trailing edge node.
  
   @param[in]   apparent_velocity   Apparent velocity.
   @param[in]   node_index          Trailing edge node index.
   
   @returns The wake emission velocity, at the node_index'th trailinge edge node.
*/
Eigen::Vector3d
LiftingSurface::wake_emission_velocity(const Eigen::Vector3d &apparent_velocity, int node_index) const{
    Vector3d wake_emission_velocity;
    
    if (Parameters::wake_emission_follow_bisector) {
        Vector3d wake_normal = wake_normals.row(node_index);
        
        // Project apparent velocity onto wake emission plane:
        wake_emission_velocity = -(apparent_velocity - apparent_velocity.dot(wake_normal) * wake_normal);
        
    } else {
        // Emit wake in direction of apparent velocity:
        wake_emission_velocity = -apparent_velocity;
        
    }
    
    // Done:
    return wake_emission_velocity;
}
