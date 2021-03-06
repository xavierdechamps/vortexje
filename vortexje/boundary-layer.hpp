//
// Vortexje -- Boundary layer base class.
//
// Copyright (C) 2014 Baayen & Heinz GmbH.
//
// Authors: Jorn Baayen <jorn.baayen@baayen-heinz.com>
//

#ifndef __BOUNDARY_LAYER_HPP__
#define __BOUNDARY_LAYER_HPP__

#include <memory>

#include <Eigen/Core>

#include <vortexje/surface.hpp>

namespace Vortexje
{

/**
   Per-surface boundary layer base class.
   
   @brief Boundary layer base class.
*/
class BoundaryLayer
{
public:
    /**
       Destructor.
    */
    virtual ~BoundaryLayer() {};
    
    /**
       Normally, this function solves the relevant boundary layer equations.  Here, it does nothing.

       @param[in]   freestream_velocity   Freestream velocity vector.
       @param[in]   surface_velocities    (n x 3)-matrix of surface velocities.
       
       @returns true on success.
    */    
    virtual bool recalculate(const Eigen::Vector3d &freestream_velocity, const Eigen::MatrixXd &surface_velocities) = 0;
    
    /**
       Returns the thickness of the boundary layer at the given panel.
       
       @param[in]   surface   Reference surface.
       @param[in]   panel     Reference panel.
   
       @returns Boundary layer thickness of the given panel.
    */
    virtual double thickness(const std::shared_ptr<Surface> &surface, int panel) const = 0; 
    
    /**
       Returns the velocity in the boundary layer at the given panel, at the given wall distance.
       
       @param[in]   surface   Reference surface.
       @param[in]   panel     Reference panel.
       @param[in]   y         Wall distance.
   
       @returns Boundary layer velocity of the given panel, at the given wall distance.
    */
    virtual Eigen::Vector3d velocity(const std::shared_ptr<Surface> &surface, int panel, double y) const = 0; 
    
    /**
       Returns the blowing velocity for the given panel.
   
       @param[in]   surface   Reference surface.
       @param[in]   panel     Reference panel.
   
       @returns Blowing velocity for the given panel.
    */
    virtual double blowing_velocity(const std::shared_ptr<Surface> &surface, int panel) const = 0;
    
    /**
       Returns the friction force acting on the given panel.
   
       @param[in]   surface   Reference surface.
       @param[in]   panel     Reference panel.
   
       @returns Friction force acting on the given panel.
    */
    virtual Eigen::Vector3d friction(const std::shared_ptr<Surface> &surface, int panel) const = 0;
};

};

#endif // __BOUNDARY_LAYER_HPP__
