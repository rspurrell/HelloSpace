 /*

  __  __          ___    ___                
 /\ \/\ \        /\_ \  /\_ \               
 \ \ \_\ \     __\//\ \ \//\ \     ___      
  \ \  _  \  /'__`\\ \ \  \ \ \   / __`\    
   \ \ \ \ \/\  __/ \_\ \_ \_\ \_/\ \L\ \   
    \ \_\ \_\ \____\/\____\/\____\ \____/   
     \/_/\/_/\/____/\/____/\/____/\/___/    
 ____                                      
/\  _`\                                    
\ \,\L\_\  _____      __      ___     __   
 \/_\__ \ /\ '__`\  /'__`\   /'___\ /'__`\ 
   /\ \L\ \ \ \L\ \/\ \L\.\_/\ \__//\  __/ 
   \ `\____\ \ ,__/\ \__/.\_\ \____\ \____\
    \/_____/\ \ \/  \/__/\/_/\/____/\/____/
             \ \_\                         
              \/_/              by Reaktor (http://hellospace.reaktor.com/)
 */
return (function (window, Controls, Quaternion, Vec3) {
  "use strict";
  const vectorPrecision = 5e-3; //1e-6
  const anglePrecision = 0.0175; // within 1 degree (was 0.001)
  const gravConst = 6.7408e-11;
  const Constants = {
    gravitation: 6.7e-11,
    timeStep: 1.0 / 60.0,
    maxThrust: 100,
    gravitationFalloff: 1.25,
    maxGimbalAngleDegrees: 45,
    fuelConsumption: 0.25,
    rocketSize: [4, 1.5, 1.5],
    rcsThrustRatio: 0.05
  };

  function YPR(y, p, r) {
    this.Yaw = y;
    this.Pitch = p;
    this.Roll = r;
  }
  
  function Utils() {};
  Utils.VecDiff = function(vec2, vec1) {
    if (!vec1) {
      vec1 = new Vec3();
    }
    return new Vec3(vec2.x - vec1.x, vec2.y - vec1.y, vec2.z - vec1.z);
  };
  Utils.VecDot = function (vec1, vec2) {
    return vec1.x * vec2.x + vec1.y * vec2.y + vec1.z * vec2.z;
  };
  Utils.VecLgth = function(vec) {
    return Math.sqrt(vec.x * vec.x + vec.y * vec.y + vec.z * vec.z);
  };
  Utils.VecLgthSqr = function(vec) {
    return vec.x * vec.x + vec.y * vec.y + vec.z * vec.z;
  };
  Utils.VecMult = function(s, vec) {
    return new Vec3(vec.x * s, vec.y * s, vec.z * s);
  };
  Utils.VecUnit = function(vec) {
    var l = Utils.VecLgth(vec);
    if (!l) {
      return new Vec3();
    }
    return new Vec3(vec.x / l, vec.y / l, vec.z / l);
  };
  Utils.ToDirectionVec = function(ra) {
    return new Vec3(
      ra.magnitude * Math.sin(ra.polar) * Math.cos(ra.azimuth), // x
      ra.magnitude * Math.sin(ra.polar) * Math.sin(ra.azimuth), // y
      ra.magnitude * Math.cos(ra.polar) // z
    );
  };
  Utils.GetApsides = function (m, v, r) {
    // see http://space.stackexchange.com/questions/1904/how-to-programmatically-calculate-orbital-elements-using-position-velocity-vecto
    var ev, evNorm, sma, E, u, vNorm = Utils.VecLgth(v), rNorm = Utils.VecLgth(r);
    u = Constants.gravitation * m;
    E = (vNorm * vNorm / 2) - (u / rNorm);
    sma = -u / 2 * E;
    //ev = 1 / mu * (((norm(v) ** 2 - mu / norm(r)) * r) - (dot(r, v) * v));
    ev = Utils.VecMult(1 / u, Utils.VecDiff(Utils.VecMult((vNorm * vNorm) - (u / rNorm), r), Utils.VecMult(Utils.VecDot(r, v), v)));
    evNorm = Utils.VecLgth(ev);
    return { apoapsis: sma * (1 + evNorm), periapsis: sma * (1 - evNorm) };
  };
  Utils.CalculateOrbitalData = function(rPos, rVel, pPos, pVel, pMass) {// see https://en.wikipedia.org/wiki/Orbital_elements
    // Relative Position Vector
    var r = Utils.VecDiff(rPos, pPos);

    // Distance between bodies
    var rNorm = Utils.VecLgth(r);

    // Relative Velocity Vector
    var v = Utils.VecDiff(rVel, pVel);

    // Relative speed between bodies
    var vNorm = Utils.VecLgth(v);

    // Specific Angular Momentum
    var h = r.cross(v);

    // Specific Angular Momentum Magnitude
    var hNorm = Utils.VecLgth(h);

    // Standard Gravitational Parameter
    var mu = Constants.gravitation * pMass;

    // Vector to Ascending Node
    // var n = new Vec3(-h.y, h.x, 0);
    var n = new Vec3(0, 0, 1).cross(h);

    // Eccentricity Vector
    var ev = Utils.VecMult(1 / mu, Utils.VecDiff(Utils.VecMult((vNorm * vNorm) - (mu / rNorm), r), Utils.VecMult(Utils.VecDot(r, v), v)));

    // Specific Orbital Energy
    var E = (vNorm * vNorm / 2) - (mu / rNorm);

    // Eccentricity
    var e = Utils.VecLgth(ev);

    // Semi-major axis
    var a = -mu / (2 * E);

    // Apoapsis
    var ap = a * (1 + e);

    // Periapsis
    var pe = a * (1 - e);

    var arg_pe, // ω: Argument of Periapsis
        f, // θ or f: True Anomaly
        i = 0, // Inclination: the angle between the angular momentum vector and its z component.
        raan, // Ω Omega: Right Ascension of Ascending Node / Longitude of Ascending Node
        SMALL_NUMBER = 1e-15;

    if (hNorm > 0) {
      i = Math.acos(h.z / hNorm);
    }

    if (Math.abs(i - 0) < SMALL_NUMBER) {
      // For non-inclined orbits, raan is undefined;
      // set to zero by convention
      raan = 0;
      if (Math.abs(e - 0) < SMALL_NUMBER) {
        // For circular orbits, place periapsis
        // at ascending node by convention
        arg_pe = 0;
      }
      else
      {
        // Argument of periapsis is the angle between
        // eccentricity vector and its x component.
        arg_pe = Math.acos(ev.x / e);
      }
    }
    else {
      // Right ascension of ascending node is the angle
      // between the node vector and its x component.
      raan = Math.acos(n.x / Utils.VecLgth(n));
      if (n.y < 0) {
        raan = 2 * Math.PI - raan;
      }

      // Argument of periapsis is angle between
      // node and eccentricity vectors.
      arg_pe = Math.acos(Utils.VecDot(n, ev) / (Utils.VecLgth(n) * e));
    }

    if (Math.abs(e - 0) < SMALL_NUMBER) {
      if (Math.abs(i - 0) < SMALL_NUMBER) {
        // True anomaly is angle between position
        // vector and its x component.
        f = Math.acos(r.x / rNorm);
        if (v.x > 0) {
          f = 2 * Math.PI - f;
        }
      }
      else {
        // True anomaly is angle between node
        // vector and position vector.
        f = Math.acos(Utils.VecDot(n, r) / (Utils.VecLgth(n) * rNorm));
        if (Utils.VecDot(n, v) > 0) {
          f = 2 * Math.PI - f;
        }
      }
    }
    else {
      if (ev.z < 0) {
        arg_pe = 2 * Math.PI - arg_pe;
      }

      // True anomaly is angle between eccentricity
      // vector and position vector.
      f = Math.acos(Utils.VecDot(ev, r) / (e * rNorm));
      if (Utils.VecDot(r, v) < 0) {
        f = 2 * Math.PI - f;
      }
    }

    return { Ap: ap, Pe: pe, ArgPe: arg_pe, Eccentricity: e, Inclination: i, RAAN: raan, SemiMajorAxis: a, TrueAnomaly: f, VecToAscNode: n };
  };
  Utils.ToRelAngle = function(vec) {
    var l = Utils.VecLgth(vec);
    var theta = Math.acos(vec.z / l);
    var phi = Math.atan2(vec.y, vec.x);
    return { polar: theta, azimuth: phi, magnitude: l };
  };
  Utils.ToEulerianAngle = function(q) {
    var ysqr = q.y * q.y;

    // yaw (z-axis rotation)
    var t3 = 2.0 * (q.w * q.z + q.x * q.y);
    var t4 = 1.0 - 2.0 * (ysqr + q.z * q.z);

    // pitch (y-axis rotation)
    var t2 = 2.0 * (q.w * q.y - q.z * q.x);
    t2 = t2 > 1.0 ? 1.0 : t2;
    t2 = t2 < -1.0 ? -1.0 : t2;

    // roll (x-axis rotation)
    var t0 = 2.0 * (q.w * q.x + q.y * q.z);
    var t1 = 1.0 - (2.0 * (q.x * q.x + ysqr));

    return new YPR(Math.atan2(t3, t4), Math.asin(t2), Math.atan2(t0, t1));
  }
  Utils.ToEuler = function(q) {
    var heading, attitude, bank;
    var x = q.x, y = q.y, z = q.z, w = q.w;

    var test = x * y + z * w;
    if (test > 0.499) { // singularity at north pole
      heading = 2 * Math.atan2(x, w);
      attitude = Math.PI / 2;
      bank = 0;
    }
    else if (test < -0.499) { // singularity at south pole
      heading = -2 * Math.atan2(x, w);
      attitude = -Math.PI / 2;
      bank = 0;
    }
    else {
      let sqx = x * x;
      let sqy = y * y;
      let sqz = z * z;
      heading = Math.atan2(2 * y * w - 2 * x * z , 1 - 2 * sqy - 2 * sqz); // Heading
      attitude = Math.asin(2 * test); // attitude
      bank = Math.atan2(2 * x * w - 2 * y * z , 1 - 2 * sqx - 2 * sqz); // bank
    }

    return new YPR(attitude, heading, bank);
  };
  Utils.ToEuler2 = function(q) {
    var heading, attitude, bank;
    var sqw = q.w * q.w;
    var sqx = q.x * q.x;
    var sqy = q.y * q.y;
    var sqz = q.z * q.z;
    var unit = sqx + sqy + sqz + sqw; // if normalised is one, otherwise is correction factor
    var test = q.x * q.y + q.z * q.w;

    if (test > 0.499 * unit) { // singularity at north pole
      heading = 2 * Math.atan2(q.x, q.w);
      attitude = Math.PI / 2;
      bank = 0;
    }
    else if (test < -0.499 * unit) { // singularity at south pole
      heading = -2 * Math.atan2(q.x, q.w);
      attitude = -Math.PI/2;
      bank = 0;
    }
    else {
      heading = Math.atan2(2 * q.y * q.w - 2 * q.x * q.z, sqx - sqy - sqz + sqw);
      attitude = Math.asin(2 * test / unit);
      bank = Math.atan2(2 * q.x * q.w - 2 * q.y * q.z, -sqx + sqy - sqz + sqw)
    }

    return new YPR(heading, attitude, bank);
  };
  Utils.PointToLocalFrame = function(position, quaternion, worldPoint, result) {
    var tmpQuat = new Quaternion();
    result = result || new Vec3();
    worldPoint.vsub(position, result);
    quaternion.conjugate(tmpQuat);
    tmpQuat.vmult(result, result);
    return result;
  };
  Utils.PointToWorldFrame = function(position, quaternion, localPoint, result) {
    result = result || new Vec3();
    quaternion.vmult(localPoint, result);
    result.vadd(position, result);
    return result;
  };
  

  function AutoPilot() {
    var self = this;

    function Init() {
      self.CurrentStage = 0;
      self.IsRunning = false;
      self.IsLogging = false;
      self.Stages = [];
      self.States = [];
      self.Tick = 0;
      self.SetControls(0, 0, 0, 0);
    }
    
    function Start(state) {
      self.IsRunning = true;
      self.IsLogging = true;
      log(self.Tick + ": Launch!");
      log(self.Tick + ": Initializing " + self.Stages[self.CurrentStage].Name);
    }
    
    function Main(state) {
      Update(state);
      while (self.Stages[self.CurrentStage].Execute(state)) {
        log(self.Tick + ": " + self.Stages[self.CurrentStage].Name + " complete.");
        self.CurrentStage++;
        log(self.Tick + ": Initializing " + self.Stages[self.CurrentStage].Name);
      }
      self.Tick++;
      return new Controls(self.Controls);
    }

    function Update(state) {
      var ps;
      state.rocket.speed = Utils.VecLgth(state.rocket.velocity);
      state.rocket.forward = state.rocket.rotation.vmult(state.rocket.position).unit();
      state.rocket.angle = Utils.ToRelAngle(state.rocket.forward);
      state.rocket.axisAngle = state.rocket.rotation.toAxisAngle();

      for (var i = 0; i < state.planetStates.length; i++) {
        ps = state.planetStates[i];
        ps.direction = Utils.VecDiff(ps.position, state.rocket.position);
        ps.orbital = Utils.CalculateOrbitalData(state.rocket.position, state.rocket.velocity, ps.position, ps.velocity, ps.mass);
        ps.directionUnit = ps.direction.unit();
        ps.distanceSurface = Utils.VecLgth(ps.direction) - ps.radius;
        ps.angle = Utils.ToRelAngle(ps.direction);
        ps.speed = Utils.VecLgth(ps.velocity);
        ps.rotReq = new Quaternion();
        ps.rotReq.setFromVectors(state.rocket.forward, ps.direction);
        ps.EulerAngle = Utils.ToEulerianAngle(ps.rotReq);
      }
      self.States.push(state);
      
      /// mult rotation by position to get facing vector?
      // log(state.rocket.forward)
      // log(Utils.ToEulerianAngle(state.rocket.rotation), 10);
    }

    function log(msg, intverval) {
      if (self.IsLogging && (!intverval || self.Tick % intverval == 0)) {
        console.log(msg);
      }
    }

    self.AddStage = function(name, operations, options) {
      self.Stages.push(new Stage(name, self, log, operations, options));
    };
    
    self.SetControls = function(t, y, p, r) {
      t = t === undefined ? self.Controls.thrust : t;
      y = y === undefined ? self.Controls.rcs.yaw : y;
      p = p === undefined ? self.Controls.rcs.pitch : p;
      r = r === undefined ? self.Controls.rcs.roll : r;
      self.Controls = { thrust: t, rcs: { pitch: p, yaw: y, roll: r }};
    };

    var execute = function(state) {
      Start(state);
      execute = Main;
      return execute(state);
    };

    self.GoToMoon = function(state) {
      return execute(state);
    };

    Init();
  }

  function Stage(name, ap, log, operations, options) {
    this.AP = ap;
    this.Controls = { t: 0, p: 0, y: 0, r: 0 };
    this.IsComplete = false;
    this.Log = log;
    this.Name = name;
    this.Options = options || { Delay: false };
    this.Ticks = [];
    this.Execute = function(state) {
      operations.call(this, state);
      this.Ticks.push({
        Tick: this.AP.Tick,
        Thrust: this.Controls.t,
        Yaw: this.Controls.y,
        Pitch: this.Controls.p,
        Roll: this.Controls.r
      });
      this.AP.SetControls(this.Controls.t, this.Controls.y, this.Controls.p, this.Controls.r);
      return this.IsComplete;
    };
    this.ApplyOptions = function() {
      this.Controls.t = this.Options.t || 0;
      this.Controls.y = this.Options.y || 0;
      this.Controls.p = this.Options.p || 0;
      this.Controls.r = this.Options.r || 0;
    };
  }

  window.AutoPilot = new AutoPilot();

  window.AutoPilot.AddStage("Burn Sub-Orbital", function (state) {
    if (this.Ticks.length == 0) {
      this.Controls.t = 1;
      this.Controls.y = 0.1075;
      this.TargetFuelVolume = 60;
    }
    if (state.rocket.fuel.volume <= this.TargetFuelVolume) {
      this.IsComplete = true;
    }
  });

  window.AutoPilot.AddStage("Arrest Spin 1", ArrestAngularVelocity, { Delay: 40, y: 1, p: 0, r: 0 });

  window.AutoPilot.AddStage("Burn to orbit", function (state) {
    if (this.Ticks.length == 0) {
      this.Controls.t = 1;
      this.TargetFuelVolume = 45;
    }
    if (state.rocket.fuel.volume <= this.TargetFuelVolume) {
      this.IsComplete = true;
    }
  });

  window.AutoPilot.AddStage("Point to Moon", PointPlanet, { DelayUntil: 1000, PlanetIndex: 1 });

  window.AutoPilot.AddStage("Arrest Spin 2", ArrestAngularVelocity, { y: 1, p: 1, r: 1 });

  window.AutoPilot.AddStage("Burn to Moon", function(state) {
    if (this.Ticks.length == 0) {
      this.Controls.t = 1;
      this.TargetFuelVolume = 20;
    }

    if (state.rocket.fuel.volume <= this.TargetFuelVolume) {
      this.IsComplete = true;
    }
  });

  window.AutoPilot.AddStage("Wait for Lunar Encounter", function() {
    if (this.AP.Tick == 2390) {
      this.IsComplete = true;
    }
  });

  window.AutoPilot.AddStage("Prep Landing Angle", function(state) {
    var azi = state.rocket.angle.azimuth;
    if (this.Ticks.length == 0) {
      this.InitAzi = azi;
      this.HalfAngle = Math.PI / 3 + 0.2;
      this.Controls.y = 1;
    }

    var aziTrav = Math.abs(azi - this.InitAzi);
    if (aziTrav >= this.HalfAngle) {
      this.IsComplete = true;
    }
  });

  window.AutoPilot.AddStage("Arrest Spin 3", ArrestAngularVelocity, { y: 1, p: 1, r: 1 });

  window.AutoPilot.AddStage("Landing", function(state) {
    if (this.Ticks.length == 0) {
      this.Controls.t = 1;
      this.TargetFuelVolume = 0;
    }

    if (state.rocket.fuel.volume <= this.TargetFuelVolume) {
      this.IsComplete = true;
    }
  });

  window.AutoPilot.AddStage("Final", function (state) {
    if (this.Ticks.length == 0) {
      this.Log(this.AP.Tick + ": Engine Shutdown");
    }
  });

  function ArrestAngularVelocity(state) {
    var x = state.rocket.angularVelocity.x,
        y = state.rocket.angularVelocity.y,
        z = state.rocket.angularVelocity.z;

    if (this.Ticks.length == 0) {
      if (!this.Options.Delay) {
        this.ApplyOptions();
      }
    }
    if (this.Options.Delay) {
      if (this.Ticks.length < this.Options.Delay) {
        return;
      }
      else if (this.Ticks.length == this.Options.Delay) {
        this.ApplyOptions();
        this.Options.Delay = false;
      }
    }

    if (state.rocket.angularVelocity.almostZero(vectorPrecision)) {
      this.IsComplete = true;
      return;
    }

    this.Controls.y = Math.abs(this.Controls.y);
    if (Math.abs(z) < 0.05) {
      this.Controls.y = this.Controls.y - 0.1;
      if (this.Controls.y < 0.1) {
        this.Controls.y = 0.1;
      }
    }
    if (z - vectorPrecision > 0) {
      this.Controls.y = -this.Controls.y;
    }
    else if (z + vectorPrecision > 0) {
      this.Controls.y = 0;
    }

    this.Controls.p = Math.abs(this.Controls.p);
    if (Math.abs(y) < 0.5) {
      this.Controls.p = this.Controls.p - 0.1;
      if (this.Controls.p < 0.1) {
        this.Controls.p = 0.1;
      }
    }
    if (y - vectorPrecision > 0) {
      this.Controls.p = -this.Controls.p;
    }
    else if (y + vectorPrecision > 0) {
      this.Controls.p = 0;
    }

    this.Controls.r = Math.abs(this.Controls.r);
    if (Math.abs(x) < 0.5) {
      this.Controls.r = this.Controls.r - 0.1;
      if (this.Controls.r < 0.1) {
        this.Controls.r = 0.1;
      }
    }
    if (x - vectorPrecision > 0) {
      this.Controls.r = -this.Controls.r;
    }
    else if (x + vectorPrecision > 0) {
      this.Controls.r = 0;
    }
  }

  function PointPlanet(state) {
    var ps = state.planetStates[this.Options.PlanetIndex];
    var y = ps.EulerAngle.Yaw,
        p = ps.EulerAngle.Pitch,
        r = ps.EulerAngle.Roll,
        yAbs = Math.abs(y),
        pAbs = Math.abs(p),
        rAbs = Math.abs(r);
    
    if (this.Ticks.length == 0) {
      if (this.Options.DelayUntil) {
        this.Options.Delay = this.Options.DelayUntil - this.AP.Tick;
      }
      if (!this.Options.Delay) {
        this.ApplyOptions();
      }
      this.HalfYawAbs = yAbs / 2;
      this.HalfPitchAbs = pAbs / 2;
      this.HalfRollAbs = rAbs / 2;
    }
    if (this.Options.Delay) {
      if (this.Ticks.length < this.Options.Delay) {
        return;
      }
      else if (this.Ticks.length == this.Options.Delay) {
        this.ApplyOptions();
        this.Options.Delay = false;
      }
    }

    if ((yAbs < this.HalfYawAbs || yAbs < anglePrecision) &&
        (pAbs < this.HalfPitchAbs || pAbs < anglePrecision) &&
        (rAbs < this.HalfRollAbs || rAbs < anglePrecision)) {
      this.IsComplete = true;
      return;
    }

    this.Controls.y = 0;
    this.Controls.p = 0;
    this.Controls.r = 0;

    if (yAbs > anglePrecision && yAbs > this.HalfYawAbs) {
      this.Controls.y = y > 0 ? 1 : -1;
    }
    if (pAbs > anglePrecision && pAbs > this.HalfPitchAbs) {
      this.Controls.p = p > 0 ? 1 : -1;
    }
    // if (Math.abs(r) > this.HalfRollAbs) {
    //   this.Controls.r = r > 0 ? 1 : -1;
    // }
  }

  function PointAwayFromPlanet(state) {
    var revVec = Utils.VecDiff(state.rocket.position, state.planetStates[this.Options.PlanetIndex].position);
    var revRot = new Quaternion();
    revRot.setFromVectors(state.rocket.forward, revVec);
    var revEulerAngle = Utils.ToEulerianAngle(revRot);
    var y = revEulerAngle.Yaw,
        p = revEulerAngle.Pitch,
        r = revEulerAngle.Roll,
        yAbs = Math.abs(y),
        pAbs = Math.abs(p),
        rAbs = Math.abs(r);
    
    if (this.Ticks.length == 0) {
      if (this.Options.DelayUntil) {
        this.Options.Delay = this.Options.DelayUntil - this.AP.Tick;
      }
      if (!this.Options.Delay) {
        this.ApplyOptions();
      }
      this.HalfYawAbs = yAbs / 2;
      this.HalfPitchAbs = pAbs / 2;
      this.HalfRollAbs = rAbs / 2;
    }
    if (this.Options.Delay) {
      if (this.Ticks.length < this.Options.Delay) {
        return;
      }
      else if (this.Ticks.length == this.Options.Delay) {
        this.ApplyOptions();
        this.Options.Delay = false;
      }
    }

    if ((yAbs < this.HalfYawAbs || yAbs < anglePrecision) &&
        (pAbs < this.HalfPitchAbs || pAbs < anglePrecision) &&
        (rAbs < this.HalfRollAbs || rAbs < anglePrecision)) {
      this.IsComplete = true;
      return;
    }

    this.Controls.y = 0;
    this.Controls.p = 0;
    this.Controls.r = 0;

    if (yAbs > anglePrecision && yAbs > this.HalfYawAbs) {
      this.Controls.y = y > 0 ? 1 : -1;
    }
    if (pAbs > anglePrecision && pAbs > this.HalfPitchAbs) {
      this.Controls.p = p > 0 ? 1 : -1;
    }
    // if (Math.abs(r) > this.HalfRollAbs) {
    //   this.Controls.r = r > 0 ? 1 : -1;
    // }
  }

  return window.AutoPilot.GoToMoon;
})(window, Controls, Quaternion, Vec3);
