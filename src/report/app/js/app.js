'use strict';

function copy_properties(source, target) {
    for (var prop in source) {
        if (source.hasOwnProperty(prop)) {
            if (source[prop] != null) {
                target[prop] = source[prop];
            }
        }
    }
}

function Dummy(obj) { copy_properties(obj, this); }
function ReportIndexEntry(obj) {
    copy_properties(obj, this);
    this.get_event = function() {
        if (this.event_best != null) {
            return this.event_best;
        } else {
            return this.event_reference;
        }
    }

    this.get_solution_type = function() {
        if (this.event_best != null) {
            return 'best';
        } else {
            return 'ref';
        }
    }
}

var yaml_type_map = [
    ['!grond.ReportIndexEntry', ReportIndexEntry],
    ['!pf.Event', Dummy],
    ['!pf.MomentTensor', Dummy],
    ['!grond.ParameterStats', Dummy],
    ['!grond.TargetBalancingAnalyserResult', Dummy],
    ['!grond.ResultStats', Dummy],
    ['!grond.WaveformMisfitTarget', Dummy],
    ['!grond.WaveformMisfitConfig', Dummy],
    ['!grond.WaveformTargetGroup', Dummy],
    ['!grond.PhaseRatioTarget', Dummy],
    ['!grond.PhaseRatioTargetGroup', Dummy],
    ['!grond.FeatureMeasure', Dummy],
    ['!grond.CMTProblem', Dummy],
    ['!pf.MTSource', Dummy],
    ['!pf.HalfSinusoidSTF', Dummy],
    ['!grond.PlotCollection', Dummy],
    ['!grond.PlotGroup', Dummy],
    ['!grond.PlotItem', Dummy],
    ['!grond.PNG', Dummy],
    ['!grond.PDF', Dummy],
];

function make_constructor(type) {
    var type = type;
    var construct = function(data) {
        return new type(data);
    };
    return construct;
}

var yaml_types = [];
for (var i=0; i<yaml_type_map.length; i++) {
    var type = yaml_type_map[i][1]
    var t = new jsyaml.Type(yaml_type_map[i][0], {
        kind: 'mapping',
        instanceOf: type,
        construct: make_constructor(type),
    });
    yaml_types.push(t);
}

var report_schema = jsyaml.Schema.create(yaml_types);

function parse_fields_float(fields, input, output, error, factor) {
    parse_fields(fields, input, output, error, factor, parseFloat);
}

function parse_fields_int(fields, input, output, error) {
    parse_fields(fields, input, output, error, 1.0, parseInt);
}

function parse_fields(fields, input, output, error, factor, parse) {
    for (var i=0; i<fields.length; i++) {
        var field = fields[i];
        if (input[field].length == 0) {
            val = null;
        } else {
            var val = parse(input[field]) * factor;
            if (val.isNaN) {
                error[field] = true;
                return false;
            }
        }
        output[field] = val;
    }
}


angular.module('reportApp', ['ngRoute'])

    .config(function($routeProvider, $locationProvider) {
        $locationProvider.hashPrefix('');
        $routeProvider
            .when('/', {
                controller: 'ReportListController',
                templateUrl: 'report_list.html',
            })
            .when('/:report_path*/', {
                controller: 'ReportController',
                templateUrl:'report.html',
            })
            .otherwise({
                template: '<div class="container">Not found!</div>',
            });
    })

    .factory('YamlDoc', function($http) {

        var funcs = {};
        funcs.query = function(path, loaded, options) {
            $http.get(path, {'responseType': 'text'}).then(
                function(response) {
                    var doc = jsyaml.safeLoad(response.data, options);
                    loaded(doc);
                }
            );
        };
        return funcs;
    })

    .factory('YamlMultiDoc', function($http) {

        var funcs = {};
        funcs.query = function(path, loaded, options) {
            $http.get(path, {'responseType': 'text'}).then(
                function(response) {
                    jsyaml.safeLoadAll(response.data, loaded, options);
                }
            );
        };
        return funcs;
    })

    .controller('NavigationController', function($scope, $route, YamlDoc, YamlMultiDoc, $routeParams, $location) {
        $scope.$route = $route;
        $scope.$location = $location;
        $scope.$routeParams = $routeParams;

        $scope.active = function(path) {
            return (path === $location.path().substr(0,path.length)) ? 'active' : '';
        };

    })

    .controller('ReportListController', function($scope, YamlMultiDoc) {
        var report_entries = [];

        var get_order_key_funcs = {
            'event': function(x) {
                return x.get_event().name;
            },
            'time': function(x) {
                return x.get_event().time;
            },
            'magnitude': function(x) {
                return x.get_event().magnitude;
            },
            'depth': function(x) {
                return x.get_event().depth;
            },
            'type': function(x) {
                return x.get_solution_type() + '.' + x.get_event().name;
            },
            'run': function(x) {
                return x.problem_name;
            }};

        var order_skey = 'event';
        var get_order_key = get_order_key_funcs[order_skey];

        $scope.set_order = function(skey) {
            order_skey = skey;
            get_order_key = get_order_key_funcs[skey];
        };

        var ordered_lines = {};

        $scope.get_ordered_report_entries = function() {
            if (order_skey in ordered_lines === false) {
                var lines = [];
                for (var i=0; i<report_entries.length; i++) {
                    var line = {
                        bg_class: null,
                        report_entry: report_entries[i]};

                    lines.push(line);
                }

                lines.sort(function(a, b) { 
                    var a_key = get_order_key(a.report_entry);
                    var b_key = get_order_key(b.report_entry);

                    if (a_key < b_key) {
                        return -1;
                    }
                    if (a_key > b_key) {
                        return 1;
                    }
                    return 0;
                });

                var i_bg_class = 0;
                var name = '';
                var bg_classes = ['even', 'odd'];
                for (var i=0; i<lines.length; i++) {
                    var rname = lines[i].report_entry.get_event().name;
                    if (rname != name) {
                        i_bg_class++;
                        name = rname;
                    }
                    lines[i].bg_class = bg_classes[i_bg_class % 2];
                }
                ordered_lines[order_skey] = lines;
            }
            return ordered_lines[order_skey];
        };

        YamlMultiDoc.query(
            'report_list.yaml',
            function(doc) { report_entries.push(doc); ordered_lines = {}; },
            {schema: report_schema});
    })

    .controller('ReportController', function(
            $scope, YamlDoc, YamlMultiDoc, $routeParams, $location, $anchorScroll) {

        $scope.stats = null;
        $scope.plot_groups = [];
        $anchorScroll.yOffset = 60;

        $scope.path = $routeParams.report_path;

        var plot_group_path = function(group_ref) {
            return $scope.path + '/plots/' + group_ref[0] + '/' + group_ref[1] + '/' + group_ref[0] + '.' + group_ref[1];
        };

        $scope.image_path = function(group, item) {
            return plot_group_path([group.name, group.variant]) + '.' + item.name + '.d100.png'
        }

        YamlDoc.query(
            $scope.path + '/stats.yaml',
            function(doc) { $scope.stats = doc; },
            {schema: report_schema});


        var query_group = function(group_ref) {
            YamlDoc.query(
                plot_group_path(group_ref) + '.plot_group.yaml',
                function(doc) { $scope.plot_groups.push(doc); },
                {schema: report_schema});
        };

        YamlMultiDoc.query(
            $scope.path + '/plots/plot_collection.yaml',
            function(doc) {
                Array.forEach(
                    doc.group_refs,
                    query_group);
                },
            {schema: report_schema});

        $scope.scrollTo = function(id) {
            var old = $location.hash();
            $location.hash(id);
            $anchorScroll();
            //reset to old to keep any additional routing logic from kicking in
            $location.hash(old);
        };
    })

    .filter('eround', function() {
        return function(input, std) {
            if (std > 0) {
                var ndig = - Math.floor(Math.log10(std)) + 1;
                var factor = Math.pow(10, ndig);
                return Math.round(input * factor) / factor;
            } else {
                return input;
            }
        };
    })

    .filter('dotalign', function() {
        return function(input) {
            input = input.toString();
            var dotpos = input.indexOf('.');
            if (dotpos == -1) {
                dotpos = input.length;
            }
            var fill = ' ';
            return fill.repeat(Math.max(0, 5 - dotpos)) + input;
        };
    })

    .run(function($rootScope, $location, $anchorScroll, $routeParams) {
      //when the route is changed scroll to the proper element.
      $rootScope.$on('$routeChangeSuccess', function(newRoute, oldRoute) {
        $location.hash($routeParams.scrollTo);
        $anchorScroll();  
      });
    });
