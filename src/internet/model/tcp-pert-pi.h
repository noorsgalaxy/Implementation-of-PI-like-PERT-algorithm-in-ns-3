/* -*- Mode:C++; c-file-style:"gnu"; indent-tabs-mode:nil; -*- */
/*
 * Copyright (c) 2016 ResiliNets, ITTC, University of Kansas
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License version 2 as
 * published by the Free Software Foundation;
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
 * Author: Truc Anh N. Nguyen <annguyen@ittc.ku.edu>
 *
 * James P.G. Sterbenz <jpgs@ittc.ku.edu>, director
 * ResiliNets Research Group  http://wiki.ittc.ku.edu/resilinets
 * Information and Telecommunication Technology Center (ITTC)
 * and Department of Electrical Engineering and Computer Science
 * The University of Kansas Lawrence, KS USA.
 */


#ifndef TCPPERTRED_H
#define TCPPERTRED_H

#include "ns3/nstime.h"
#include "ns3/boolean.h"
#include "ns3/data-rate.h"
#include "ns3/timer.h"
#include "ns3/event-id.h"
#include "ns3/random-variable-stream.h"
#include "ns3/tcp-congestion-ops.h"
#include "ns3/tcp-recovery-ops.h"

namespace ns3 {
////  Delete the below vegas Related comment And add tcpPertPi related comments


/**
 * \ingroup congestionOps
 *
 * \brief An implementation of TCP Vegas
 *
 * TCP Vegas is a pure delay-based congestion control algorithm implementing a proactive
 * scheme that tries to prevent packet drops by maintaining a small backlog at the
 * bottleneck queue.
 *
 * Vegas continuously measures the actual throughput a connection achieves as shown in
 * Equation (1) and compares it with the expected throughput calculated in Equation (2).
 * The difference between these 2 sending rates in Equation (3) reflects the amount of
 * extra packets being queued at the bottleneck.
 *
 *              actual = cwnd / RTT        (1)
 *              expected = cwnd / BaseRTT  (2)
 *              diff = expected - actual   (3)
 *
 * To avoid congestion, Vegas linearly increases/decreases its congestion window to ensure
 * the diff value fall between the 2 predefined thresholds, alpha and beta.
 * diff and another threshold, gamma, are used to determine when Vegas should change from
 * its slow-start mode to linear increase/decrease mode.
 *
 * Following the implementation of Vegas in Linux, we use 2, 4, and 1 as the default values
 * of alpha, beta, and gamma, respectively.
 *
 * More information: http://dx.doi.org/10.1109/49.464716
 */

class TcpPertPi : public TcpNewReno
{
public:
  /**
   * \brief Get the type ID.
   * \return the object TypeId
   */
  static TypeId GetTypeId (void);

  /**
   * Create an unbound tcp socket.
   */
  TcpPertPi (void);

  /**
   * \brief Copy constructor
   * \param sock the object to copy
   */
  TcpPertPi (const TcpPertPi& sock);
  virtual ~TcpPertPi (void);

  virtual std::string GetName () const;

  /**
   * \brief Compute RTTs needed to execute Vegas algorithm
   *
   * The function filters RTT samples from the last RTT to find
   * the current smallest propagation delay + queueing delay (minRtt).
   * We take the minimum to avoid the effects of delayed ACKs.
   *
   * The function also min-filters all RTT measurements seen to find the
   * propagation delay (baseRtt).
   *
   * \param tcb internal congestion state
   * \param segmentsAcked count of segments ACKed
   * \param rtt last RTT
   *
   */
  virtual void PktsAcked (Ptr<TcpSocketState> tcb, uint32_t segmentsAcked,
                          const Time& rtt);

  /**
   * \brief Enable/disable Vegas algorithm depending on the congestion state
   *
   * We only start a Vegas cycle when we are in normal congestion state (CA_OPEN state).
   *
   * \param tcb internal congestion state
   * \param newState new congestion state to which the TCP is going to switch
   */
  virtual void CongestionStateSet (Ptr<TcpSocketState> tcb,
                                   const TcpSocketState::TcpCongState_t newState);

  /**
   * \brief Adjust cwnd following Vegas linear increase/decrease algorithm
   *
   * \param tcb internal congestion state
   * \param segmentsAcked count of segments ACKed
   */
  virtual void IncreaseWindow (Ptr<TcpSocketState> tcb, uint32_t segmentsAcked);

  /**
   * \brief Get slow start threshold following Vegas principle
   *
   * \param tcb internal congestion state
   * \param bytesInFlight bytes in flight
   *
   * \return the slow start threshold value
   */
  virtual uint32_t GetSsThresh (Ptr<const TcpSocketState> tcb,
                                uint32_t bytesInFlight);

  virtual Ptr<TcpCongestionOps> Fork ();


private:

  void CheckChangeLossProb(Ptr<TcpSocketState> tcb);
  void CalculateP();
  void CheckAndSetAlpha(Ptr<TcpSocketState> tcb);
  void UpdatePertVars(const Time& l_rtt);


  Time m_thresh1;                   //!<
  Time m_thresh2;                   //!<
  Time m_thresh3;                   //!<
  double m_maxp;                  //!<probability value at thresh2_
  double m_beta;                  //!< Decresing Factor
  double m_alpha;                 //!<Maximum alpha value set to 32
  double m_dProb;                 //!<weighted average of drop probability (p) computed in TFRC style
  double m_erProb;                //!<weighted average of early response probability (p') computed in TFRC style
  double m_pertProb;                //!<
  double m_nd;                    //!<number of packets sent in the current since last packet loss
  double m_ndp1;                  //!<number of packets sent since last early response (for the estimation of p')
  Time m_minRtt;  
  Time m_maxRtt;

  double m_pertSrtt;  ////double
  uint32_t m_changeWindow;
  uint32_t m_competeRegionCounter;
  uint32_t  m_highspeedRegionCounter;
  uint32_t  m_mode;
  double m_maxProb; 
     
  Time  m_lastUpdatedAlpha;      ////check this it will be int or in double // Time
  double m_alphaMax;
  bool m_sender;
  double m_this;
  double(m_a);
double(m_b);
double(m_qRef);
double(m_qOld);
  EventId m_rtrsEvent;
  SequenceNumber32 m_lastEarlyResponseSeq;
  std::vector<double> m_weight;
  std::vector<double> m_historyND;
  std::vector<double> m_historyNDp1;
   


};

} // namespace ns3

#endif // TCPPERTRED_H
