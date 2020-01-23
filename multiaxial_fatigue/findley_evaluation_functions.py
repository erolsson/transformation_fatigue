# Import standard python modules
import numpy as np
import time, sys, pickle, multiprocessing
from math import sin, cos, pi, sqrt


def evaluate_findley(combined_stress, a_cp, worker_run_out_time=3600, chunk_size=100,
                     cpus=multiprocessing.cpu_count(),
                     w_pool=None, search_grid=5, **_):
    s_time = time.time()
    if not w_pool:
        worker_pool = multiprocessing.Pool(processes=cpus)
    else:
        worker_pool = w_pool

    # Get size of the read data
    load_steps, rows, columns = combined_stress.shape
    print " Read %2i load steps with %i stress tensors" % (load_steps, rows)

    # Create workload vector
    work_loads = range(0, rows, chunk_size)
    if not work_loads[-1] == rows:
        work_loads.append(rows)
    print " Number or work pieces to process: ", len(work_loads)

    # Create storage point for Findley results
    fatigue_results = np.empty(work_loads[-1], dtype=float)
    if True:
        # Submit workloads for evaluation of Findley stress
        print " Computing critical plane stress:"
        findley_load_step_jobs = []
        for work_load in range(0, len(work_loads) - 1):
            findley_load_step_jobs.append(worker_pool.apply_async(findley_worker,
                                                                  [[a_cp[work_loads[work_load]:work_loads[work_load + 1]],
                                                                    combined_stress[:, work_loads[work_load]:
                                                                                    work_loads[work_load + 1], :],
                                                                    search_grid]]))
        # Retrieve results for workloads
        for work_load, findley_load_step_job in enumerate(findley_load_step_jobs):
            # print "Working with workload " + str(work_load)

            fatigue_results[work_loads[work_load]:work_loads[work_load + 1]] = findley_load_step_job.get(
                worker_run_out_time)
            print ".",
            # print "Done with workload " + str(work_load)
            sys.stdout.flush()  # Force output of buffered content
        print "\n Done, Total Time: %1.2f" % (time.time() - s_time)

        # # Save to pickle dump file    
        # print " Pickle Findly stresses to file,",
        # pickle_handle = open('Findly2.pkl', 'wb')
        # pickle.dump(fatigue_results.shape,pickle_handle)
        # pickle.dump(fatigue_results,pickle_handle)
        # pickle_handle.close()
        # print "Done!"
    if not w_pool:
        worker_pool.close()
        worker_pool.join()
    return fatigue_results


# ----------------------------------------------------------------------------------------------------------------------


def eval_findley(a_cp, stress_matrix, search_grid, mod=False):
    def smallest_enclosing_circle(xp, yp, xo=None, yo=None):
        # [xc, yc, R] = smallest_enclosing_circle(X, Y)
        #
        # Purpose:
        # Calculate the Smallest Enclosing Circle of a set of points
        #
        # Input:  
        # xp                       np array with X coordinates of the points, size: 1 x nr_points
        # yp                       np array with Y coordinates of the points, size: 1 x nr_points
        #
        # Input used for recursive use only:
        # xo                  array with 1, 2 or 3 X coordinates of outermost points, size: 1 x nr_outer_points
        # yo                  array with 1, 2 or 3 Y coordinates of outermost points, size: 1 x nr_outer_points
        #
        # Output:
        # xc                 X coordinate of smallest enclosing circle
        # yc                 Y coordinate of smallest enclosing circle
        # radius             radius of smallest enclosing circle
        #
        # History: 
        # 14-Dec-2006 creation by FSta
        #             based on an example by Yazan Ahed (yash78@gmail.com),
        #             who based his code on a Java applet by Shripad Thite (http://heyoka.cs.uiuc.edu/~thite/mincircle/)

        # Initialize xo, yo, nr_outer_points
        if xo is None or yo is None:
            xo = []
            yo = []
            nr_outer_points = 0
        else:
            nr_outer_points = np.size(xo)

        # Compute new center point coordinates and radius
        if nr_outer_points == 0:
            xc = np.mean(xp)
            yc = np.mean(yp)
            radius = 0
        elif nr_outer_points == 1:
            xc = xo[0]
            yc = yo[0]
            radius = 0
        elif nr_outer_points == 2:
            xc = (xo[0] + xo[1]) / 2
            yc = (yo[0] + yo[1]) / 2
            radius = sqrt((xo[0] - xc) ** 2 + (yo[0] - yc) ** 2)

        elif nr_outer_points == 3:
            xc = (xo[2] ** 2 * (yo[0] - yo[1]) + (xo[0] ** 2 + (yo[0] - yo[1]) * (yo[0] - yo[2])) *
                  (yo[1] - yo[2]) + xo[1] ** 2 * (-yo[0] + yo[2])) / (2 * (xo[2] * (yo[0] - yo[1]) +
                                                                           xo[0] * (yo[1] - yo[2]) + xo[1] *
                                                                           (-yo[0] + yo[2])))
            yc = (yo[1] + yo[2]) / 2 - (xo[2] - xo[1]) / (yo[2] - yo[1]) * (xc - (xo[1] + xo[2]) / 2)
            radius = sqrt((xo[0] - xc) ** 2 + (yo[0] - yc) ** 2)
            return xc, yc, radius
        else:
            print "warning... caught an unexpected mode... in elif 1"
            raise RuntimeError

        # Check if points are within the circle
        for i, _ in enumerate(xp):
            if (xp[i] - xc) ** 2 + (yp[i] - yc) ** 2 > radius ** 2:
                if xo is not None or yo is not None:
                    if not xp[i] in xo or not yp[i] in yo:
                        if nr_outer_points == 0:
                            xo = [xp[i]]
                            yo = [yp[i]]
                        elif nr_outer_points == 1:
                            xo = [xo[0], xp[i]]
                            yo = [yo[0], yp[i]]
                        elif nr_outer_points == 2:
                            xo = [xo[0], xo[1], xp[i]]
                            yo = [yo[0], yo[1], yp[i]]
                        else:
                            print "warning... caught an unexpected mode... in elif 2"
                            raise RuntimeError

                        [xc, yc, radius] = smallest_enclosing_circle(xp[0:i + 1], yp[0:i + 1], xo, yo)

        return xc, yc, radius

    def get_transform_matrix(theta_deg, phi_deg):
        # Radians
        theta_r = pi * theta_deg / 180.0
        phi_r = pi * phi_deg / 180.0

        # Multiaxial fatigue, Marquis, Eq 1.3 & 1.5
        a11 = cos(theta_r) * sin(phi_r)
        a12 = sin(theta_r) * sin(phi_r)
        a13 = cos(phi_r)
        a21 = -sin(theta_r)
        a22 = cos(theta_r)
        a23 = 0
        a31 = -cos(theta_r) * cos(phi_r)
        a32 = -sin(theta_r) * cos(phi_r)
        a33 = sin(phi_r)

        # Compose transformation matrix
        trans_matrix = np.array([[a11 ** 2, a12 ** 2, a13 ** 2, 2 * a11 * a12, 2 * a11 * a13, 2 * a13 * a12],
                                 [a21 ** 2, a22 ** 2, a23 ** 2, 2 * a21 * a22, 2 * a21 * a23, 2 * a23 * a22],
                                 [a31 ** 2, a32 ** 2, a33 ** 2, 2 * a31 * a32, 2 * a31 * a33, 2 * a33 * a32],
                                 [a11 * a21, a12 * a22, a13 * a23, a11 * a22 + a12 * a21, a13 * a21 + a11 * a23,
                                  a12 * a23 + a13 * a22],
                                 [a11 * a31, a12 * a32, a13 * a33, a11 * a32 + a12 * a31, a13 * a31 + a11 * a33,
                                  a13 * a32 + a12 * a33],
                                 [a21 * a31, a22 * a32, a23 * a33, a21 * a32 + a22 * a31, a23 * a31 + a21 * a33,
                                  a22 * a33 + a23 * a32]])
        return trans_matrix

        #     Search Space

    phi_space = 90
    theta_space = 180

    #     Get shape of stress matrix
    loadsteps, points, no_stress_components = stress_matrix.shape

    # Temporary vector for storage of derived nodal stresses
    # nodal_temp_array=np.empty([loadsteps,3], dtype=float)

    # Result array [theta, phi, max_sigma_n, max_tau_amplitude, F]
    findley_vec = np.zeros((points, 5))

    # Loop over all planes 
    first_run = True  # The first time, do not compare just store data
    for theta in range(0, theta_space + search_grid, search_grid):
        for phi in range(-phi_space, phi_space + search_grid, search_grid):

            # Compute the transformation matrix for the considered plane
            q = get_transform_matrix(theta, phi)

            j = 0  # Iterator
            # For the currently considered planed, evaluate sigma_n, tau_1, tau_2, 
            # Delta_tau, F for every node for the load history (time domain)
            for node_s_hist_vector in np.rollaxis(stress_matrix, 1):  # Loop over the stress history for a specific node
                # Compute shear stress and normal stresses for the load history
                s_prim = np.dot(node_s_hist_vector, q.T)
                # Evaluate the smallest enclosing circle to get the shear stress amplitude 
                # (i.e. the radius of the circle)
                x, y, max_tau_amplitude = smallest_enclosing_circle(s_prim[:, 3], s_prim[:, 4])

                # Evaluate the largest normal stress on the plane for the load history
                max_sigma_n = s_prim[:, 0].max()

                # Store result if the are larger than the current value (except for first plane, then just store data)
                if first_run:
                    findley_vec[j, :] = theta, phi, max_sigma_n, max_tau_amplitude, (max_tau_amplitude +
                                                                                     a_cp[j] * max_sigma_n)
                else:
                    if findley_vec[j, 4] < max_tau_amplitude + a_cp[j] * max_sigma_n:
                        findley_vec[j, 4] = max_tau_amplitude + a_cp[j] * max_sigma_n
                        findley_vec[j, 0] = phi
                        findley_vec[j, 1] = theta
                        findley_vec[j, 2] = max_sigma_n
                        findley_vec[j, 3] = max_tau_amplitude
                    # Correct implementation above by erolsson
                    # if findley_vec[j,2] < max_sigma_n:
                    #    findley_vec[j,2] = max_sigma_n
                    # if findley_vec[j,3] < max_tau_amplitude:
                    #    findley_vec[j,0] = phi
                    #    findley_vec[j,1] = theta
                    #    findley_vec[j,3] = max_tau_amplitude

                j += 1

            if first_run:
                first_run = False

                # Evaluate Findley stress
    #    F = max_sigma_n  +  a_cp * max_tau_amplitude
    # findley_vec[:,4]=findley_vec[:,2]+a_cp*findley_vec[:,3]

    return findley_vec[:, 4]  # Just return findley stress


# ----------------------------------------------------------------------------------------------------------------------


def findley_worker(job_arguments):
    # Expand recieved arguments
    # a_cp,stress_matrix=job_arguments

    # Create results array

    try:
        return eval_findley(job_arguments[0], job_arguments[1], job_arguments[2])

    except KeyboardInterrupt:
        raise KeyboardInterrupt
    except:
        print "Problem with computing Findley stress:"
        print sys.exc_info()[0]  # - Exit type:', sys.exc_info()[0]
        print sys.exc_info()[1]  # - Exit type:', sys.exc_info()[0]
        raise
